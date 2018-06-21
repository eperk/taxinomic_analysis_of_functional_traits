r_fossilrix <- as.matrix(royer_data_fossil_int)
r_fossilrix <- r_fossilrix[,1:3]

locate.fossil(royer_tree, r_fossilrix)

#####temptesting
royer_data_temp <- royer_data_sub[-c(3,5)]
temp_fossil <- all_fossil[1,]
temp_full <- rbind(royer_data_temp, temp_fossil)
temp_full <- temp_full[-c(1)]

temp_full <- as.matrix(temp_full)

fossil_mean <- aggregate(all_fossil[,-1], by=list(all_fossil$binomial), mean)

temptree <- locate.fossil(royer_tree, temp_full, search="exhaustive",plot=TRUE)

add.arrow(tree=temptree,tip="Fagopsis_longifolia",col="red",lwd=3,hedl=0.06,angle=50)


######testingcontall_fossil$binomial
######located.fossil w/ date constraints and edge constrants?
###### should just work
###### ask revell for help
#### get_mrca
##### figtree
tree_quercus <- add.species.to.genus(tree_plant, "Quercus_sp", where=c("root"))
tree_royer_q <- add.species.to.genus(royer_tree, "Quercus_sp", where = c("root"))
tree_fossil <- add.species.to.genus(tree_plant, fossil_mean$Group.1, where=c("root"))



plot(tree_quercus, type = "fan", show.tip.label = FALSE)
plot(tree_royer_q, type = "fan")
add.arrow(tree=tree_royer_q,tip="Quercus_sp",col="blue")

