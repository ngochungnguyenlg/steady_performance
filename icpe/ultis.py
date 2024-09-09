import numpy as np
from icpe.config import *

def estimate_warmup_time(fork_m, w, wi):
        wt = 0
        sit = 0
        sic = 0
        for e in fork_m:
            ni = int(e/100)
            it = ni*e
            sit += it
            if sit >w:
                wt += sit
                sit = 0
                sic += 1
            if sic >= wi:
                break
        return wt

def select_performance_measurements(wt, r, i, M):
    sit = 0
    sic = 0
    t = 0
    Mconf = []
    for e in M:
        ni = e / 100 
        it = ni * e
        if t >= wt:
            sit += it
            Mconf.append(e)
            if sit > r:
                sit = 0
                sic += 1
            if sic >= i:
                break
        t += e
    return Mconf


def filter_outliers(data, window_size=200):
    """
    Lọc các outlier từ dữ liệu sử dụng cửa sổ trượt (sliding window) 
    và điều kiện median ± 3 × (99%ile − 1%ile).
    
    Args:
    - data (list hoặc np.array): Dữ liệu đầu vào.
    - window_size (int): Kích thước cửa sổ trượt.

    Returns:
    - np.array: Dữ liệu đã được lọc bỏ các outlier.
    """
    data = np.array(data)
    window_size = min(window_size, len(data))  # Đảm bảo window_size không lớn hơn độ dài của data
    noops = np.ceil(0.1 / np.mean(data))
    skip = int(np.ceil(window_size / np.log1p(noops)))

    # Khởi tạo danh sách các chỉ số ban đầu, bỏ qua những phần đầu theo skip
    indexes = list(range(skip, len(data)))
    outliers = set()

    def sliding_window(seq, window_size):
        """Hàm hỗ trợ để tạo ra các cửa sổ trượt từ danh sách các chỉ số."""
        for i in range(len(seq) - window_size + 1):
            yield seq[i:i + window_size]

    def create_isoutlier(series, window):
        """Tạo hàm để xác định xem một giá trị có phải là outlier hay không."""
        series_ = series[window]
        median = np.median(series_)
        p1, p99 = [np.quantile(series_, q) for q in [0.01, 0.99]]
        threshold = 3 * (p99 - p1)

        def _isoutlier(i):
            value = series[i]
            return abs(value - median) > threshold

        return _isoutlier

    # Kiểm tra từng cửa sổ và xác định outliers
    for window in sliding_window(indexes, window_size):
        isoutlier = create_isoutlier(data, window)
        outliers.update([i for i in window if isoutlier(i)])

    # Lọc bỏ các outlier từ dữ liệu
    filtered_data = np.setdiff1d(np.arange(len(data)), list(outliers))
    
    return data[filtered_data]


# def filter_outliers(data, window_size=200):
#     """
#     using median ± 3 × (99%ile − 1%ile)
#     Args:
#     - data (list or np.array): input data.
#     - window_size (int): size of a slicing windown

#     Returns:
#     - np.array: filter data.
#     """
#     data = np.array(data)
#     filtered_data = []
    
#     for i in range(len(data)):
#         start = max(0, i - window_size // 2)
#         end = min(len(data), i + window_size // 2 + 1)
        
#         window = data[start:end]
        
#         median = np.median(window)
#         p99 = np.percentile(window, 99)
#         p1 = np.percentile(window, 1)
#         threshold = 3 * (p99 - p1)
        
#         if abs(data[i] - median) <= threshold:
#             filtered_data.append(data[i])
    
#     return np.array(filtered_data)


import ruptures as rpt

def calculate_segs_cost(tau, penalty, model):
    tau = [0] + tau
    total_cost = 0
    
    for i in range(1, len(tau)):
        start = tau[i-1]+1
        end = tau[i]        
        if start < end-1:
            try:
                cost = model.cost.error(start, end-1) + penalty
                total_cost += cost 
            except Exception as e:
                print("Error happen, because of {}, start at {}, end at {}".format(e, start, end-1))
    return total_cost

def crops_algorithm(data, beta_min, beta_max, method = "l2", jump=5, min_size=2, max_iterator = 1000):
    '''
    params:
    
    '''
    if len(data) == 0:
        return {}
    model = rpt.Pelt(model=method, min_size=min_size, jump=jump).fit(data)
    def run_pelt(penalty):
        return model.predict(pen=penalty)    
            
    beta_set = {(beta_min, beta_max)}
    try:
        segmentations = {beta_min: run_pelt(beta_min), beta_max: run_pelt(beta_max)}
    except Exception as e:
        print (e, data)
    n= max_iterator
    while beta_set and n>0:
        # idx = np.random.randint(0, len(beta_set))
        beta0, beta1 = beta_set.pop()
        try:
            segs_0 = run_pelt(beta0)
            segs_1 = run_pelt(beta1)
            m0 = len(segs_0)
            m1 = len(segs_1)
            if m0 >m1 + 1:
                beta_int = (calculate_segs_cost(segs_1,beta1, model) - calculate_segs_cost(segs_0,beta0, model)) /(m0-m1)
                
                segments_int = run_pelt(beta_int)
                
                if len(segments_int) != m1:
                    beta_set.update([(beta0, beta_int), (beta_int, beta1)])
                    segmentations[beta_int] = segments_int
                print(len(beta_set), beta_int, m0, m1, beta0, beta1)
        except Exception as e:
            print (e, data)
        n -=1
    return segmentations


def find_elbow(penalties, costs):
    from kneed import KneeLocator
    kneedle = KneeLocator(penalties, costs, curve='convex', direction='decreasing')
    elbow_point = kneedle.elbow
    return elbow_point


def calculate_mean_and_ci(data):
    """Calculate the mean and 95% confidence interval of the mean."""
    mean = np.mean(data)
    ci = 1.96 * np.std(data, ddof=1) / np.sqrt(len(data))
    return mean, ci

def is_within_tolerance(last_mean, last_ci, new_mean, new_ci, tolerance=0.05):
    """Check if mean2 is within the tolerance of mean1."""
    
    last_low = last_mean-last_ci
    last_up = last_mean+last_ci
    
    new_low = new_mean-new_ci
    new_up = new_mean+new_ci
    
    low_difer = abs(last_low-new_low)/last_low
    up_difer = abs(last_up-new_up)/last_up
    
    return low_difer<tolerance and up_difer < tolerance

def find_stable_segment(data, segments):
    """Find stable segments based on the 5% tolerance on the CI of the mean."""
    stable_segments = []
    unstable_segments = []
    if len(segments) < 2:
        return stable_segments
    if config["using_last_stable_segmet"]==False:
        max_length_seg = -float("inf")
        max_length_seg_idx = -1
        #ignore the first segnment
        for idx in range(1,len(segments)-1):
            if segments[idx] - segments[idx+1] > max_length_seg:
                max_length_seg = segments[idx] - segments[idx+1]
                max_length_seg_idx = idx
        mean, last_ci = calculate_mean_and_ci(data[segments[max_length_seg_idx]:segments[max_length_seg_idx+1]])     
    else:
        mean, last_ci = calculate_mean_and_ci(data[segments[-2]:segments[-1]])
    
    start_stabel_range = 0
    end_stabel_range = 0
    for i in range(len(segments) - 2, -1, -1):
        segment_data = data[segments[i]:segments[i+1]]
        last_mean, ci = calculate_mean_and_ci(segment_data)
        
        if is_within_tolerance(mean, ci, last_mean,last_ci):
            stable_segments.append((segments[i], segments[i+1]))
            start_stabel_range = segments[i]
            if end_stabel_range < segments[i+1]:
                end_stabel_range = segments[i+1]
        else:
            unstable_segments.append((segments[i], segments[i+1]))
    is_consistent = True
    if start_stabel_range != end_stabel_range:
        for seg in unstable_segments:
            if seg[0] < end_stabel_range and seg[0] > start_stabel_range:
                is_consistent = False
                
    return stable_segments, is_consistent

if __name__ == "__main__":
    np.random.seed(0)
    N = np.random.randint(100,200,5)
    noise = np.random.uniform(-3,3,5)
    data = np.concatenate([np.random.normal(loc+noise[loc], 0.5, size=N[loc]) for loc in range(5)], axis=0)
    data = data.flatten()
    # Chạy thuật toán CROPS
    beta_min = 4
    beta_max = 10^5
    segmentations = crops_algorithm(data, beta_min, beta_max, max_iterator=10, method="l1")

    import matplotlib.pyplot as plt
    plt.plot(data)
    # In kết quả
    for beta, segments in segmentations.items():
        for cp in segments:
            plt.axvline(x=cp, color='r', linestyle='--')
            
    plt.show()