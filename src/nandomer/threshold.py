import statistics



def find_threshold(values, threshold = "moderate"):
    
    median = statistics.median(values)
    mad = (statistics.median([abs(x - median) for x in values]))*1.4826
    
    if threshold == "moderate":
        return median - mad
    
    if threshold == "high":
        return median
    
    if threshold == "superhigh":
        return median + mad
    raise ValueError("Invalid mode. Choose 'moderate', 'high', or 'superhigh'.")

def determine_threshold(values, threshold = "moderate"):
    if not values:
        return 0.0
    return find_threshold(values, threshold)
