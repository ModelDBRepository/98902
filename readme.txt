This is the model associated with the paper

Anderson WS, Kudela P, Cho J, Bergey GK, Franaszczuk PJ (2007) Studies
of stimulus parameters for seizure disruption using neural network
simulations. Biol Cybern 97:173-94

To use first type make in the HopkinsModel directory.

Currently I'm (Stan Anderson) running it on a 32-bit 16-node cluster,
also running an older Redhat version. I've been compiling it on a
64-bit machine, but in 32-bit mode, the 64-bit machine is running a
newer Suse version.

As far as the individual c-codes:

mknode.c is the code that sets up the intracolumnar connections and
produces a synapse table that is passed to each node.

mklink.c is the code that sets up the extracolumnar connections, which
produces an unsorted link table that undergoes further
manipulations. This code also calculates which neurons in the
simulation will undergo action potentials during a stimulation pulse.

sublink.c assigns the previously computed extracolumnar links to
specific neurons, and strips away intranodal links, turning that
information into synaptic structures as created in
mknode.c. Extranodal links (i.e. computer to computer) or maintained
as link structures.

link_sort.c sorts the link structures into an order more easily read
by the simulation code.

crl.c actually writes the link information into a file of link
structures that is passed to each node.

netclustwacnmdadiff.c is the main simulation code running on each node
with the integration loop, and accounting of synaptic currents.

clust_cn.c - handles communication via nodes via packet protocols
(very slow in this setting)

the lnet.h and clust_cn.h files are the most important header files,
and most variables can be found here.

mkpar.c enables the user to fill the parameter tables for each cell
class in terms of individual channel kinetic properties.
