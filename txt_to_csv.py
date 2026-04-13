import pandas as pd
import click

@click.command()
@click.option('--input_txt', help='The path to the text file')
@click.option('--output_csv', help='The path to the output CSV file')

def txt_to_csv(input_txt, output_csv):
    # Define the column names
    column_names = ['k', 'K(k)', 'Kcsr(k)', 'K/Kcsr(k)', 'RND(k)']

    # Read the text file into a DataFrame with the specified column names
    df = pd.read_csv(input_txt, delim_whitespace=True, comment='#', names=column_names)

    # Save the DataFrame as a CSV file
    df.to_csv(output_csv, index=False)

if __name__ == '__main__':
    txt_to_csv()