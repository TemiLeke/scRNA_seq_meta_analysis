import synapseclient
import synapseutils

syn = synapseclient.Synapse()
syn.login('tadeoye_usf','V_DKSvsA3VHDFsL')
fastq_files = synapseutils.syncFromSynapse(syn, 'syn32164389',  path='../data/mathys_pfc/fastq/')
#count_files = synapseutils.syncFromSynapse(syn, 'syn18681734', path='../data/mathys_pfc/fastq/')