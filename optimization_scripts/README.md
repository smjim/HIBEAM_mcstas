## How to run venetian blinds generator
1. Determine input codes and constant variables and input into simulation manually
2. Determine input variables and run code with those as arguments:
example:
`python3 optimization_scripts/focused_Venbla_xy_to_off.py 10 0.5 0.0005 4 55 [outDir] --detpos -0.3 -0.1`
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

