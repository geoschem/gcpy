#!/usr/bin/env python3
"""
Prints information about a GEOS-Chem Classic horizontal grid
"""
from gcpy.grid import make_grid_ll


def print_vals(values):
    """
    Prints a list of numbers 8 columns wide with 5 decimal places.

    Args
    values : list-like : Values to be displayed
    """
    cols = 8
    for i in range(0, len(values), cols):
        line = "".join(f"{x:11.5f}" for x in values[i:i+cols])
        print(line)


def create_grid_and_print_info(grid_params):
    """
    Creates a grid object and prints its metadata.

    Args
    grid_params : tuple : Resolution, lon range, lat range
    """

    # Arguments
    res = grid_params[0]
    lon_range = grid_params[1]
    lat_range = grid_params[2]
    grid_type = grid_params[3]

    # Make resolution a pretty string
    res_display = res.replace("x", "\u00B0 x ") + "\u00B0"

    # Create a grid object
    out_extent = [lon_range[0], lon_range[1], lat_range[0], lat_range[1]]
    grid = make_grid_ll(res, out_extent=out_extent)

    # Print metadata
    print(f"Name            : {res_display} {grid_type}")
    print(f"Resolution      : {res}")
    print(f"Longitude Range : {lon_range}")
    print(f"Latitude Range  : {lat_range}")
    print("Longitude centers")
    print_vals(grid["lon"])
    print("Latitude centers")
    print_vals(grid["lat"])
    print("Longitude edges")
    print_vals(grid["lon_b"])
    print("Latitude edges")
    print_vals(grid["lat_b"])


def get_grid_parameters(selection):
    """
    Returns metadata about a given grid

    Args
    selection   : int   : Number of the selected grid

    Returns
    grid_params : tuple : Resolution, lon range, lat range

    """
    #    Resolution       # Lon range            # Lat range
    grids = [
        ("4.0x5.0",       [-180.0,    180.0   ], [-90.0,  90.0 ], "global"),
        ("2.0x2.5",       [-180.0,    180.0   ], [-90.0,  90.0 ], "global"),
        ("0.5x0.625",     [-180.0,    180.0   ], [-90.0,  90.0 ], "global"),
        ("0.5x0.625",     [  60.0,    150.0   ], [-11.0,  55.0 ], "nested AS"),
        ("0.5x0.625",     [ -30.0,     50.0   ], [ 30.0,  70.0 ], "nested EU"),
        ("0.5x0.625",     [-140.0,    -40.0   ], [ 10.0,  70.0 ], "nested NA"),
        ("0.25x0.3125",   [-180.0,    180.0   ], [-90.0,  90.0 ], "global"),
        ("0.25x0.3125",   [ -20.0,     52.8125], [-37.0,  40.0 ], "nested AF"),
        ("0.25x0.3125",   [  70.0,    140.0   ], [ 32.75, 61.25], "nested AS"),
        ("0.25x0.3125",   [ -15.0,     40.0   ], [ 15.0,  55.0 ], "nested EU"),
        ("0.25x0.3125",   [ -20.0,     70.0   ], [ 12.0,  44.0 ], "nested ME"),
        ("0.25x0.3125",   [-130.0,    -60.0   ], [  9.75, 60.0 ], "nested NA"),
        ("0.25x0.3125",   [ 110.0,    180.0   ], [-50.0,   5.0 ], "nested OC"),
        ("0.25x0.3125",   [ -87.8125, -31.25  ], [-59.0,  16.0 ], "nested SA"),
        ("0.25x0.3125",   [  20.0,    180.0   ], [ 41.0,  83.0 ], "nested RU"),
        ("0.125x0.15625", [-180.0,    180.0   ], [-90.0,  90.0 ], "global"),
        ("0.125x0.15625", [ -20.0,     52.8125], [-59.0,  16.0 ], "nested AF"),
        ("0.125x0.15625", [  70.0,    140.0   ], [ 15.0,  55.0 ], "nested AS"),
        ("0.125x0.15625", [ -15.0,     40.0   ], [ 32.75, 61.25], "nested EU"),
        ("0.125x0.15625", [-130.0,    -60.0   ], [  9.75, 60.0 ], "nested NA"),
        ("0.125x0.15625", [ -87.8125, -31.25  ], [-59.0,  16.0 ], "nested SA"),
    ]

    return grids[selection]


def display_menu():
    """
    Displays a menu of grids to choose from.

    Returns
    selection : int : Number of the grid selected by the user
    """
    msg = """\
Please select a GEOS-Chem Classic horizontal grid:

1.  4.0\u00B0   x 5.0\u00B0     global

2.  2.0\u00B0   x 2.5\u00B0     global

3.  0.5\u00B0   x 0.625\u00B0   global
4.  0.5\u00B0   x 0.625\u00B0   nested AS
5.  0.5\u00B0   x 0.625\u00B0   nested EU
6.  0.5\u00B0   x 0.625\u00B0   nested NA

7   0.25\u00B0  x 0.3125\u00B0  global
8.  0.25\u00B0  x 0.3125\u00B0  nested AF
9.  0.25\u00B0  x 0.3125\u00B0  nested AS
10. 0.25\u00B0  x 0.3125\u00B0  nested EU
11. 0.25\u00B0  x 0.3125\u00B0  nested ME
12. 0.25\u00B0  x 0.3125\u00B0  nested NA
13  0.25\u00B0  x 0.3125\u00B0  nested OC
14. 0.25\u00B0  x 0.3125\u00B0  nested SA
15. 0.25\u00B0  x 0.3125\u00B0  nested RU

16. 0.125\u00B0 x 0.15625\u00B0 global
17. 0.125\u00B0 x 0.15625\u00B0 nested AF
18. 0.125\u00B0 x 0.15625\u00B0 nested AS
19. 0.125\u00B0 x 0.15625\u00B0 nested EU
20. 0.125\u00B0 x 0.15625\u00B0 nested NA
21. 0.125\u00B0 x 0.15625\u00B0 nested SA
    """
    print(msg)

    while True:
        try:
            selection = int(input("Enter a selection >>> "))
            if 1 <= selection <= 21:
                break   # valid, exit loop
            print("❌ Entry out of range. Try again.")
        except ValueError:
            print("❌ Not a number. Try again.")

    return selection - 1


def main():
    """
    Main program
    """

    # Display a list of grids and get the user's selection
    selection = display_menu()

    # Get the parameters that define the grid
    grid_params = get_grid_parameters(selection)

    # Create the grid object and print metadata
    create_grid_and_print_info(grid_params)


if __name__ == '__main__':
    main()
