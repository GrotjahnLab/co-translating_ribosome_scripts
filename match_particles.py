import starfile
import pandas as pd
import click

@click.command()
@click.option('--original_csv', type=click.Path(exists=True), help='Path to the original CSV file')
@click.option('--filtered_starfile', type=click.Path(exists=True), help='Path to the filtered star file')
@click.option('--output_csv', default='output.csv', help='Name of the output CSV file')

def match_particles(original_csv, filtered_starfile, output_csv):
    #load the original csv file
    original_df = pd.read_csv(original_csv)
    # Read the filtered star file
    star_filter = starfile.read(filtered_starfile)
    #use rlnImageName in star_filter to select particles in original_df
    selected_df = original_df[original_df['rlnImageName'].isin(star_filter['rlnImageName'])]
    #save the selected particles to a new csv file
    selected_df.to_csv(output_csv, index=False)

if __name__ == '__main__':
    match_particles()
