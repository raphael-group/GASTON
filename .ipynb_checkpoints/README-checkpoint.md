# SpatialNN

loading olfactory model:

nhS=50 # USE THIS - best looking layers
partition='lower'
optimizer='adam'
trial=3

folder='/n/fs/ragr-research/projects/network-mutations/manifold-alignment/olfactory_glmpca/intermediate_NN_v2/'
folder+=f'nhS_{nhS}_optimizer_{optimizer}_partition_{partition}_trial_{trial}/'

mod=torch.load(folder+'model_epoch_20000.pt')

loading cerebellum model:

nhS=200
trial=0
optimizer='adam'
partition='all'

folder='/n/fs/ragr-research/projects/network-mutations/manifold-alignment/slideseq_cerebellum/intermediate_NN_v2/'
folder+=f'nhS_{nhS}_optimizer_{optimizer}_partition_{partition}_trial_{trial}/'


mod=torch.load(folder+'model_epoch_20000.pt')