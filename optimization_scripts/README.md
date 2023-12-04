## How to run venetian blinds generator
1. Determine input codes and constant variables and input into simulation manually
2. If running with individually focused blades, generate optimal source positions
	a. Run `optimization_scripts/focused_Venbla_xy_to_off.py` with `show_blade_x_scan(VB_pos, hx, hy, VB_length, no_VB_outDir, image_dir, dx=0.2)` to generate blade backpropagation scan (dx gives resolution, blade generation does linear interpolation between these datapoints)
	b. Run `optimization_scripts/display_ridge.py` on the output blade backpropagation histograms to generate optimal source position calculation for given characteristic (recommended weighted average of distributions)
	c. Run `optimization_scripts/focused_Venbla_xy_to_off.py` with command line argument of optimal source pos datafile
example:  
`python3 display_ridge.py backprop_image/histogram_x_scan_* --save opt_source_pos.dat`  
the command above generates
2. Determine input variables and run code with those as arguments:
example:  
`python3 optimization_scripts/focused_Venbla_xy_to_off.py 10 0.5 0.0005 4 55 tmp/ --detpos -0.3 -0.1 --noShow --source_pos_interpolate opt_source_pos.dat`
3. Run plotter to show all plot files in the output directory
example:  
`python3 optimization_scripts/display_pdf.py [outDir] [output].pdf`

```
usage: focused_Venbla_xy_to_off.py [-h] [--detdim detwidth detheight] [--detpos detx dety] VB_z VB_length VB_thickness VB_m det_z output_dir

Create Venetian Blinds geometry with given parameters

positional arguments:
  VB_z                  Z distance between source and center of vb
  VB_length             Z length of blades
  VB_thickness          Thickness of blades
  VB_m                  VB reflectivity m value, default=4
  det_z                 Z distance between center of vb and target
  output_dir            Output directory for calculations

options:
  -h, --help            show this help message and exit
  --detdim detwidth detheight
                        width and height of detector
  --detpos detx dety    displacement of detector position
```


### Additional examples 
Run VB configuration for given parameter configuration
`python3 optimization_scripts/focused_Venbla_xy_to_off.py 10 0.5 0.0001 4 55 tmpa/ --detpos -0.3 -0.1`
