#this whole situation is hacky

# run as sage generate_ggb.sage <sobj>, generates ggb in location of sobj

load('TL2ggb.sage')
filename = sys.argv[1]
generate_ggb(filename, db=False)