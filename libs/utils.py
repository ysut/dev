from pathlib import Path
import pandas as pd

############ Load and output ############
def load_file(input_file: str, type: str) -> pd.DataFrame:
    if type == 'csv':
        df = pd.read_csv(input_file, header=0, dtype=str)
    elif type == 'tsv':
        df = pd.read_csv(input_file, header=0, dtype=str, sep='\t')
    
    return df


def configure_output(input_file: str) -> str:
    input_path = Path(input_file) 
    output_base = f'{input_path.parent}/{input_path.stem}'

    return output_base

   
def output_tsv(df: pd.DataFrame, output: str) -> None:
    df.to_csv(output, sep='\t', index=False)

