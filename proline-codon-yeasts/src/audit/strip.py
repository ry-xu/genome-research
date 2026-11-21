import sys, re
t = open(sys.argv[1]).read()
t = re.sub(r'<[^>]+>', ' ', t); t = re.sub(r'\s+', ' ', t)
open(sys.argv[2], 'w').write(t)
