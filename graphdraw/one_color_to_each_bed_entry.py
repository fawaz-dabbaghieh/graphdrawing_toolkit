import os
import sys
import random
import pdb


color_list=['snow','ghostwhite','GhostWhite','whitesmoke','WhiteSmoke','gainsboro','floralwhite','FloralWhite','oldlace','OldLace','linen','antiquewhite','AntiqueWhite','papayawhip','PapayaWhip','blanchedalmond','BlanchedAlmond','bisque','peachpuff','PeachPuff','navajowhite','NavajoWhite','moccasin','cornsilk','ivory','lemonchiffon','LemonChiffon','seashell','honeydew','mintcream','MintCream','azure','aliceblue','AliceBlue','lavender','lavenderblush','LavenderBlush','mistyrose','MistyRose','white','black','darkslategray','DarkSlateGray','darkslategrey','DarkSlateGrey','dimgray','DimGray','dimgrey','DimGrey','slategray','SlateGray','slategrey','SlateGrey','lightslategray','LightSlateGray','lightslategrey','LightSlateGrey','gray','grey','lightgrey','LightGrey','lightgray','LightGray','midnightblue','MidnightBlue','navy','navyblue','NavyBlue','cornflowerblue','CornflowerBlue','darkslateblue','DarkSlateBlue','slateblue','SlateBlue','mediumslateblue','MediumSlateBlue','lightslateblue','LightSlateBlue','mediumblue','MediumBlue','royalblue','RoyalBlue','blue','dodgerblue','DodgerBlue','deepskyblue','DeepSkyBlue','skyblue','SkyBlue','lightskyblue','LightSkyBlue','steelblue','SteelBlue','lightsteelblue','LightSteelBlue','lightblue','LightBlue','powderblue','PowderBlue','paleturquoise','PaleTurquoise','darkturquoise','DarkTurquoise','mediumturquoise','MediumTurquoise','turquoise','cyan','lightcyan','LightCyan','cadetblue','CadetBlue','mediumaquamarine','MediumAquamarine','aquamarine','darkgreen','DarkGreen','darkolivegreen','DarkOliveGreen','darkseagreen','DarkSeaGreen','seagreen','SeaGreen','mediumseagreen','MediumSeaGreen','lightseagreen','LightSeaGreen','palegreen','PaleGreen','springgreen','SpringGreen','lawngreen','LawnGreen','green','chartreuse','mediumspringgreen','MediumSpringGreen','greenyellow','GreenYellow','limegreen','LimeGreen','yellowgreen','YellowGreen','forestgreen','ForestGreen','olivedrab','OliveDrab','darkkhaki','DarkKhaki','khaki','palegoldenrod','PaleGoldenrod','lightgoldenrodyellow','LightGoldenrodYellow','lightyellow','LightYellow','yellow','gold','lightgoldenrod','LightGoldenrod','goldenrod','darkgoldenrod','DarkGoldenrod','rosybrown','RosyBrown','indianred','IndianRed','saddlebrown','SaddleBrown','sienna','peru','burlywood','beige','wheat','sandybrown','SandyBrown','tan','chocolate','firebrick','brown','darksalmon','DarkSalmon','salmon','lightsalmon','LightSalmon','orange','darkorange','DarkOrange','coral','lightcoral','LightCoral','tomato','orangered','OrangeRed','red','hotpink','HotPink','deeppink','DeepPink','pink','lightpink','LightPink','palevioletred','PaleVioletRed','maroon','mediumvioletred','MediumVioletRed','violetred','VioletRed','magenta','violet','plum','orchid','mediumorchid','MediumOrchid','darkorchid','DarkOrchid','darkviolet','DarkViolet','blueviolet','BlueViolet','purple','mediumpurple','MediumPurple','thistle','snow1','snow2','snow3','snow4','seashell1','seashell2','seashell3','seashell4','AntiqueWhite1','AntiqueWhite2','AntiqueWhite3','AntiqueWhite4','bisque1','bisque2','bisque3','bisque4','PeachPuff1','PeachPuff2','PeachPuff3','PeachPuff4','NavajoWhite1','NavajoWhite2','NavajoWhite3','LemonChiffon1','LemonChiffon2','LemonChiffon3','LemonChiffon4','cornsilk1','cornsilk2','cornsilk3','cornsilk4','ivory1','ivory2','ivory3','ivory4','honeydew1','honeydew2','honeydew3','honeydew4','LavenderBlush1','LavenderBlush2','LavenderBlush3','LavenderBlush4','MistyRose1','MistyRose2','MistyRose3','MistyRose4','azure1','SlateBlue1','SlateBlue2','SlateBlue3','SlateBlue4','RoyalBlue1','RoyalBlue2','RoyalBlue3','RoyalBlue4','blue1','blue2','blue3','blue4','DodgerBlue1','DodgerBlue2','DodgerBlue3','DodgerBlue4','SteelBlue1','SteelBlue2','SteelBlue3','SteelBlue4','DeepSkyBlue1','DeepSkyBlue2','DeepSkyBlue3','DeepSkyBlue4','SkyBlue1','SkyBlue2','SkyBlue3','SkyBlue4','LightSkyBlue1','LightSkyBlue2','LightSkyBlue3','LightSkyBlue4','SlateGray1','SlateGray2','SlateGray3','SlateGray4','LightSteelBlue1','LightSteelBlue2','LightSteelBlue3','LightSteelBlue4','LightBlue1','LightBlue2','LightBlue3','LightBlue4','LightCyan1','LightCyan2','LightCyan3','LightCyan4','PaleTurquoise1','PaleTurquoise2','PaleTurquoise3','PaleTurquoise4','CadetBlue1','CadetBlue2','CadetBlue3','CadetBlue4','turquoise1','turquoise2','turquoise3','turquoise4','cyan1','cyan2','cyan3','cyan4','DarkSlateGray1','DarkSlateGray2','DarkSlateGray3','DarkSlateGray4','aquamarine1','aquamarine2','aquamarine3','aquamarine4','DarkSeaGreen1','DarkSeaGreen2','DarkSeaGreen3','DarkSeaGreen4','SeaGreen1','SeaGreen2','SeaGreen3','PaleGreen1','PaleGreen2','PaleGreen3','SpringGreen1','SpringGreen2','SpringGreen3','DarkOliveGreen1','DarkOliveGreen2','khaki1','khaki2','khaki3','LightGoldenrod1','LightGoldenrod2','LightGoldenrod3','LightYellow1','LightYellow2','LightYellow3','LightYellow4','RosyBrown1','RosyBrown2','RosyBrown3','RosyBrown4','IndianRed1','burlywood1','burlywood2','burlywood3','wheat1','wheat2','wheat3','wheat4','LightSalmon1','LightSalmon2','DebianRed','DeepPink1','DeepPink2','DeepPink3','HotPink1','HotPink2','HotPink3','HotPink4','pink1','pink2','pink3','pink4','LightPink1','LightPink2','LightPink3','LightPink4','PaleVioletRed1','PaleVioletRed2','PaleVioletRed3','maroon1','maroon2','maroon3','VioletRed1','VioletRed2','VioletRed3','magenta1','magenta2','magenta3','magenta4','orchid1','orchid2','orchid3','orchid4','plum1','plum2','plum3','plum4','MediumOrchid1','MediumOrchid2','MediumOrchid3','MediumOrchid4','DarkOrchid1','DarkOrchid2','DarkOrchid3','DarkOrchid4','purple1','purple2','purple3','purple4','MediumPurple1','MediumPurple2','MediumPurple3','MediumPurple4','thistle1','thistle2','thistle3','thistle4','gray100','darkgrey','DarkGrey','darkgray','DarkGray','darkblue','DarkBlue','darkcyan','DarkCyan','darkmagenta','DarkMagenta','darkred','DarkRed','lightgreen','LightGreen']



if len(sys.argv) < 2:
	print("You need to give the csv file")
	sys.exit()

in_csv = sys.argv[1]
# the problem is that the CSV file is not ordered by bed entry number
# so it could have entry 1 then 5 then 1 again, I just want to have one color
# for all the alignments in one entry
# I can do one pass through the file
colors = set()
bed_entries = dict()
with open(in_csv, "r") as infile:
	next(infile)  # skipping the header
	for l in infile:
		l = l.strip().split(",")
		seq = l[-1]
		seq = seq.split("_")
		bed_ent = "_".join(seq[1:4])
		if bed_ent not in bed_entries:
			c = random.choice(color_list)
			while True:
				# pdb.set_trace()
				if c not in colors:
					colors.add(c)
					bed_entries[bed_ent] = c
					break
				else:
					c = random.choice(color_list)

counter = 0
with open(in_csv, "r") as infile:
	for l in infile:
		if counter == 0:
			print(l.strip())
			counter += 1
		else:
			l = l.strip().split(",")
			l[1] = bed_entries["_".join(l[2].split("_")[1:4])]
			print(",".join(l))
