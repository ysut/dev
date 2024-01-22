# Decorator

def print_filtering_count(func):
    def _wrapper(*args, **kwargs):
        print(f'Start {func.__name__}')
        pre = len(args[0])
        result = func(*args, **kwargs)
        post = len(result)
        print(f'Filtering : {pre} --> {post}\n')
        return result
    return _wrapper

