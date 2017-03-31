

Assignment Specification
-------------------------------------------------------------------------------
Avaiable points:

	Required:

	1. Connect three.js and a perlin library to render perlin heightmap terrain
	2. Be able to either (a) move around the terrain or (b) edit the terrain

	Optional:

	Dual contouring (15%)
	Concavity support by storing more than simple height-maps (15%)
	Terrain brush (10%)
	  + that allows concavity, not just height maps (10%, requires "concave")
	  + with multiple brushes (around 5% per brush)
	  + with texture brushes (10-50% depending on complexity)
	Procedural textures based on elevation (10%)
	  + and slope (10%)
	  + and wetness/peakiness/other non-local terrain features (10%)
	Automatic terrain modification
	  + make it drain (10%)
	  + hydraulic erosion (10%)
	  + spheroidal weathering (5%)
	  + scree (10%; requires erosion or weathering)
	  + ridge amplification/reduction (15%)
	3D noise
	  - as floating lumps (10%), or
	  - as cave digging (20%)
	Navigation
	  + camera/terrain collision prevention (10%)
	      +10% if caves can be entered but walls not penetrated
	  + walking+jumping with gravity (15%)
	  
	Completed:

		Dual Contouring (see note in sources)
		Concavity support
		Terrain brush - Multiple materials, 2 shapes, add/remove
		3D Noise - cave digging
-----------------------------------------------------------------------------------

## Usage:

Open up index.html on some web browser - tested with Google Chrome

A 64x64x64 chunk of terrain will be generated

2D Perlin noise generates a heightmap with a grass layer on top, then a dirt layer, then a rock layer

3D Simplex noise is used to dig cave-like shapes into the terrain

## Camera Controls:

Hold down the left mouse button and move the mouse to orbit the terrain.

Hold the left mouse button and the S key and move the move to zoom in and out

Hold down the right mouse button and move the mouse to pan

Z - toggles wireframe view

## Brush Controls:

The translucent sphere that you see can be used to edit the terrain

T - toggle brush shape between sphere and box

G - Move brush in -X direction

H - Move brush in +X direction

Y - Move brush in +Y direction

B - Move brush in -Y direction

In sphere mode:
+ to increase sphere radius
- to decrease sphere radius

In box mode:
+ scale box in +Y direction
- scale box in -Y direction
I - scale box in +Z direction
M - scale box in -Z direction
J - scale box in -X direction
K - scale box in +X direction

##Enter - performs selected operation

By default, the operation is delete (i.e. cut the brush shape out of the terrain)

The number keys change the operation type

1 - delete
2 - create dirt
3 - create rock
4 - create grass


##Sources: 

perlin.js - https://github.com/josephg/noisejs

TrackballControls.js - https://threejs.org/examples/misc_controls_trackball.html

Dual contouring concepts:
http://www.frankpetterson.com/publications/dualcontour/dualcontour.pdf
https://people.eecs.berkeley.edu/~jrs/meshpapers/SchaeferWarren2.pdf

Some of the math involved in Dual Contouring (solving the quadratic error function) is beyond my current understanding. In order to more quickly get to the point of experimenting with CSG operations and multiple materials, I decided to make use of a reference implementation.
I translated some code from this blog: http://ngildea.blogspot.com/2014/11/implementing-dual-contouring.html in order to implement Dual Contouring
I did modify much of the code to match the description given in the papers and to support multiple materials. I am willing to take reduced credit on the Dual Contouring part of the assignment since I did not come up with a lot of the code.

I used some concepts described in this blog to implement CSG operations:
https://upvoid.com/devblog/2013/07/terrain-engine-part-2-volume-generation-and-the-csg-tree/

