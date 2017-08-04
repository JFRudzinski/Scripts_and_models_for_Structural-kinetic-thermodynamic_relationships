Nproc=16
clustering="/data/isilon/rudzinski/soft_backup/Clustering-0.12_static/clustering"

ln -s dtraj_pca5D.dat coords && ${clustering} density -f coords -R 0.025 0.05 0.1 0.2 0.3 -p pop -d fe -n ${Nproc} -v && ln -s fe_0.100000 fe && ${clustering} density -f coords -r 0.1 -D fe -b nn -n ${Nproc} -v && ${clustering} density -f coords -r 0.1 -D fe -B nn -T 0. 0.1 -o clust -n ${Nproc} -v && ${clustering} network -p 500 --step 0.2 -v && ${clustering} density -f coords -i network_end_node_traj.dat -r 0.1 -D fe -B nn -o microstates -n ${Nproc} -v 

# && clustering mpp -i microstates -D fe -l 7 -v --concat-nframes 50001
