# exr2png
Convert to PNG/TGA image from OpenEXR image.
This also generate the blurred cubemap using same algorithm of AMD cubemapgen.
The alpha component in PNG has the reciprocal of the color strength.
To restore the color, divide RGB by alpha.

usage: exr2png.exe [-s scale] [-f type] [-c type] [-a angle] [-m level] [(+x|-x|+y|-y|+z|-z) left top right bottom dir] [infile] [outfile]

-s scale

    The low dynamic range scale.
    0.5 maps the range of 0-0.5 to 0-255, 2.0 maps 0-2.0 to 0-255.
    Note that higher scale reduce the maximum brightness.
    If not passed this option, the default is 1.0.

-f type

    The output format type.
    'png' or 'tga' is accepted.

-c type

    Enable cubemap mode.
    You can set 6 region with the order of +x -x +y -y +z -z.
    And output 6 files that have the postfix of '_[pn][xyz]'
    'type' is following type:
    filtered  : generate the filtered cubemap.
    irradiance: generate the irradiance cubemap.
    none      : generate the non filtered cubemap.

-a angle

    A half of the filter corn angle. the default is 1.0.
    In general, when the mip level goes up one, it will double.
 
-m level

    A maximum mipmap level. level number of texture is  generated.
    An image is reduced to 1/2 for each level. The default is 1, and the maximum is 16.

+x left top right bottom dir

-x left top right bottom dir

+y left top right bottom dir

-y left top right bottom dir

+z left top right bottom dir

-z left top right bottom dir

    Set the source cubemap region.
    The reagion is the rectangle of pixels.
    'left' and 'top', or 'right' and 'bottom' defines a half open range with 0 origin.
    'dir' indicates the top edge of an image, it may take a one of following characters:
    
      l Left is actual top edge.
      t Top is actual top edge.
      r Right is actual top edge.
      b Bottom is actual top edge.
    
    Therefore this options will have following.
    +x 0 0 256 256
    It should set all of six region. the layout is following:
    
    -------------> U direction
    |    +--+
    |    |+y|
    | +--+--+--+--+
    | |-x|+z|+x|-z|
    | +--+--+--+--+
    |    |-y|
    |    +--+
    |
    v  V direction

infile

    OpenEXR image file.

outfile

    The image file that converted from infile.
    If not passed this option, use the infile that has replaced extension to '.png.' or '.tga'.
    In cubemap mode, add some suffix, so the actual output filename is:
    
      outfile[1-16]_(px|nx|py|ny|pz|nz).(png|tga)
    
    It will generating the six image files for each level.
