#Update the parents populations, so we have an updated population of individuals as parents to have improvement

# replace 20 oldest parents with 10 new parents from EYT, and 10 parents from AYT, my example
Parents <- c(Parents[21:nParents], AYT[1:10], EYT)
# quiza esto lo pueda poner como funcion, no pq estricta necesidad, pero pa ser mas chulo
