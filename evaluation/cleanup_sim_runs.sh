cd sim_data
for i in {1..70} 
do
 rm -r d$i
done
cd ../sim_runs
for i in {1..70} 
do
 rm -r d$i
 rm run\_$i\_monitor.txt
done


