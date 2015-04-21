require(devtools)
require(roxygen2)
document()
install()
check()

stop('airsea loaded')

require(lineprof)
source('ncp_prof.R')
load('test.rdata')
test$NCPx1 = test$NCPx0 * test$timePeriod
test$com = O2NCP(test, entrainment = F, asVolume = T)
test$sim = O2NCP.simple(test, entrainment = F, asVolume = T)
test$diff = test$sim - test$com
test$diffper = test$diff / ((test$sim + test$com)/2) * 100
lhs$com = O2NCP(lhs, entrainment = F, asVolume = T)
lhs$sim = O2NCP.simple(lhs, entrainment = F, asVolume = T)
x = lineprof(O2NCP(test[1,]))
system.time(O2NCP(test))
# shine(x)

