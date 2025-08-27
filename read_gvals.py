import pickle,os,sys
from pprint import pprint

icg_dir='ICG/'

fname=icg_dir+'gvals_IH.pkl'
with open(fname, 'rb') as f:
	gvals_IH = pickle.load(f, encoding='latin1')

pprint(gvals_IH)

sys.exit()

# now let's print some stats of how many keys each gvals file has
files = [ x for x in os.listdir(icg_dir) if x.startswith('gvals') ]
for fname in files:
	with open(icg_dir+fname, 'rb') as f:
		gvals = pickle.load(f, encoding='latin1')
	print(fname, len(gvals.keys()))
