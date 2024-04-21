import pandas as pd

## Decorator
def print_filtering_count(func):
    def _wrapper(*args, **kwargs):
        print(f'Start {func.__name__}')
        pre = len(args[0])
        result = func(*args, **kwargs)
        post = len(result)
        print(f'Filtering : {pre} --> {post}\n')
        return result

    return _wrapper


## Filtering functions
@print_filtering_count
def extract_snv(df: pd.DataFrame, key: str) -> pd.DataFrame:
    df = df[df[key] == 'snv']

    return df


@print_filtering_count
def extract_denovo(df: pd.DataFrame) -> pd.DataFrame:
    df = df.astype({'vqslod': float, 'denovogear': float,
                    'triodenovo': float, 'dnmfilter': float})
    df = df[(df['vqslod'] > -7.18)
            & ((df['denovogear'] > 0.02)
            | (df['triodenovo'] > 5.72)
            | (df['dnmfilter'] > 0.196)
            | (df['denovofilter'] == 'TRUE'))]

    return df


@print_filtering_count
def exclude_intergenic(df: pd.DataFrame, key: str) -> pd.DataFrame:
    df = df[df[key] != 'intergenic_region']

    return df


@print_filtering_count
def exclude_utr(df: pd.DataFrame, key: str) -> pd.DataFrame:
    df = df[(df[key] != '5_prime_UTR_variant')
            & (df[key] != '5_prime_UTR_premature_start_codon_gain_variant')
            & (df[key] != '3_prime_UTR_variant')]

    return df


@print_filtering_count
def exclude_no_transcripts(df: pd.DataFrame, key: str) -> pd.DataFrame:
    df = df[df[key] != 'intragenic_variant']

    return df


@print_filtering_count
def exclude_up_down_stream(df: pd.DataFrame, key: str) -> pd.DataFrame:
    df = df[(df[key] != 'upstream_gene_variant')
            & (df[key] != 'downstream_gene_variant')]

    return df


@print_filtering_count
def exclude_tf_binding(df: pd.DataFrame, key: str) -> pd.DataFrame:
    df = df[df[key] != 'TF_binding_site_variant']

    return df


@print_filtering_count
def exclude_truncating_var(df: pd.DataFrame, key: str) -> pd.DataFrame:
    df = df[(df[key] != 'stop_gained')
            & (df[key] != 'stop_lost&splice_region_variant')]
    
    return df


@print_filtering_count
def exclude_start_lost(df: pd.DataFrame, key: str) -> pd.DataFrame:
    df = df[(df[key] != 'start_lost')
            & (df[key] != 'start_lost_splice_region_variant')]

    return df


def all_exclude(df: pd.DataFrame) -> pd.DataFrame:
    extract_snv(df)
    exclude_intergenic(df)
    exclude_utr(df)
    exclude_no_transcripts(df)
    exclude_up_down_stream(df)
    exclude_tf_binding(df)

    return df
