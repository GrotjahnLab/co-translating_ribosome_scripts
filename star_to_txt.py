import starfile
import pandas as pd
import click

@click.command()
@click.option('--input_star', help='The path to the star file')
@click.option('--output_txt', help='The path to the output txt file')
def star_to_txt(input_star, output_txt):
    # Read the star file
    star = starfile.read(input_star)

    # Save the columns named as rlnCoordinateX, rlnCoordinateY, rlnCoordinateZ into a new data frame
    coordinates = star[['rlnCoordinateX', 'rlnCoordinateY', 'rlnCoordinateZ']]

    # Save the data frame into a .txt file without the title of the columns
    coordinates.to_csv(output_txt, index=False, header=False, sep='\t')

if __name__ == '__main__':
    star_to_txt()

