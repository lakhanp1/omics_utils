
library(data.table)



######################################################################
## Example 1

## these are very good examples of reshaping the data frame
df1 <- data.frame(
  id = 1:10,
  time = as.Date('2009-01-01') + 0:9,
  Q3.2.1 = rnorm(10, 0, 1),
  Q3.2.2 = rnorm(10, 0, 1),
  Q3.2.3 = rnorm(10, 0, 1),
  Q3.3.1 = rnorm(10, 0, 1),
  Q3.3.2 = rnorm(10, 0, 1),
  Q3.3.3 = rnorm(10, 0, 1)
)

## solution1: using reshape()
## slower
ans1 = reshape(data = df1, direction = "long", 
               varying = list(c("Q3.2.1", "Q3.2.2", "Q3.2.3"), c("Q3.3.1", "Q3.3.2", "Q3.3.3")),
               v.names = c("Q3.2", "Q3.3"),
               times = c("Q1", "Q2", "Q3")
)

## solution2: using data.table::melt()
## faster
ans2 = melt.data.table(as.data.table(df1), id.vars = c("id", "time"), 
                       measure.vars = list(c("Q3.2.1", "Q3.2.2", "Q3.2.3"), c("Q3.3.1", "Q3.3.2", "Q3.3.3")), 
                       value.name = c("Q3.2", "Q3.3"))


## benchmarking
microbenchmark(
  ans1 = reshape(data = df1, direction = "long", 
                 varying = list(c("Q3.2.1", "Q3.2.2", "Q3.2.3"), c("Q3.3.1", "Q3.3.2", "Q3.3.3")),
                 v.names = c("Q3.2", "Q3.3"),
                 times = c("Q1", "Q2", "Q3")
  ),
  ans2 = melt.data.table(as.data.table(df1), id.vars = c("id", "time"), 
                         measure.vars = list(c("Q3.2.1", "Q3.2.2", "Q3.2.3"), c("Q3.3.1", "Q3.3.2", "Q3.3.3")), 
                         value.name = c("Q3.2", "Q3.3"))
)




######################################################################
## Example 2
df2 = data.frame(student=c(rep(1,5),rep(2,5)), month=c(1:5,1:5),  
                 quiz1p1=seq(20,20.9,0.1),quiz1p2=seq(30,30.9,0.1),  
                 quiz2p1=seq(80,80.9,0.1),quiz2p2=seq(90,90.9,0.1))      


## solution1: 
ans1 = reshape(df2, direction="long", idvar=c("student", "month"),
               varying = list(c("quiz1p1", "quiz2p1"), 
                              c("quiz1p2", "quiz2p2")), 
               v.names = c("p1", "p2"), times = c("quiz1", "quiz2"))



