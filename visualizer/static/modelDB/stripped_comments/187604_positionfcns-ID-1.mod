NEURON {
  SUFFIX nothing
}


FUNCTION get_x_pos(gid, gmin, BinNumX, BinNumYZ, binSizeX) {
	LOCAL CellNum, tmp
	CellNum=gid - gmin+1
	tmp = floor((CellNum-1)/BinNumYZ)
	get_x_pos =  fmod(tmp,BinNumX)*binSizeX+binSizeX/2.0
	
}

FUNCTION get_y_pos(gid, gmin, BinNumY, BinNumZ, binSizeY) {
	LOCAL CellNum, tmp, pos
	CellNum=gid - gmin+1
	tmp = floor((CellNum-1)/BinNumZ)
	get_y_pos =  fmod(tmp,BinNumY)*binSizeY+binSizeY/2.0
}

FUNCTION get_z_pos(gid, gmin, BinNumZ, binSizeZ, ZHeight) {
	LOCAL CellNum, pos
	CellNum=gid - gmin+1
	get_z_pos = fmod((CellNum-1),BinNumZ)*binSizeZ+binSizeZ/2+ZHeight
}