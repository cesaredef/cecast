## Cesare de Filippo
## Some functions for cecast (CoEstimate human Contamination And Split Time)

###################################

# Sequenctial matching for string.
# A cheap version which take the best match after comparing each charachter in target with the equivalent (i.e. at the same position) in options. 
seq_match <- function(target, options) {
	target_s <- tolower(strsplit(target, split="")[[1]])
	options_s <- lapply(strsplit(options, split=""), tolower)
	matches <- sapply(options_s, function(z) sum(z[1] == target_s[1]))
	for(i in seq_along(target_s)[-1]) {
		matches = matches + sapply(options_s, function(z) sum(z[i] == target_s[i], na.rm=T))
	}
	if(sum(matches) == 0) {
		return(NA) # no match
	} else {
		best_match <- options[which.max(matches)]
		return(best_match)
	}
}

###################################
## Run the simulations with scrm ##

## It can be improved but it works.
sim_with_scrm = function(nsims, basepairs=1e6, RHO=52, THETA=58, BASEPAIRS=1e6, rescale=TRUE,
	split_time=NULL, split_pop="v", N_test=2, N_vindija=2, N_chagyrskaya=2, N_altai=2, N_denisova=2, N_europe=2, N_africa=2, pops = c("v", "c", "a", "d", "x", "y", "o"),
	run=TRUE, params=NULL, cmd_append=NULL) {
	# The order of the populations are vindija (1 | v), chagyrskaya (2 | c), altai (3 | a), denisova (4 | d), europe (5 | y), africa (6 | x) and chimp (7 | o). Notice that chimp is not present in the output.
	test_pop = rep(0, length(pops))
	test_pop[which(pops == split_pop)] = N_test
	if(is.null(split_time)) {
		stop("split_time must be specified.")
	}
	test_pop = paste("-eI", split_time, paste(test_pop, collapse=" "))
	## Rescale the values RHO and THETA depending on the ratio basepairs/BASAPAIRS
	# Keep in mind that the simulatiuons are scaled with:
	# mu <- 1.45e-8; N0 <- 1000; basepairs <- 1e6
	# so THETA is 1.45e-8*4*1000*1e6
	if (rescale) {
		RHO = RHO*basepairs/BASEPAIRS
		THETA = THETA*basepairs/BASEPAIRS
	}
	# Fixed parameters from Chagyrskaya paper.
	# Also t_x = 300kya, a whatever fixed paparmeters
	par_to_sim = list(t_HND=4.793103, t_ND=3.655172, t_Nd=1.176724, t_Nc=0.787931,  t_XY=0.948276, s_Nv=0.446552, s_Nc=0.676724, s_Nd=1.051724, s_D=0.620690, tm_NH=0.60345, m_NH=0.975)
	if (!is.null(params)) {
		if( sum(names(params) %in% names(par_to_sim)) == length(params) ) {
			par_to_sim[names(params)] = params
		} else {
			stop("all names of params should match.")
		}
	}
	N_tot = sum(as.numeric(c(N_test, N_vindija, N_chagyrskaya, N_altai, N_denisova, N_europe, N_africa)))
	main = sub("N_europe", N_europe, sub("N_africa", N_africa, paste(N_tot, nsims, "-t", THETA, "-r", RHO, basepairs, "-I 7 0 0 0 0 N_europe N_africa 0 -l 100r -SC abs ")))
	main1 = "-eI s_Nv N_vindija 0 0 0 0 0 0 \
        -eI s_Nc 0 N_chagyrskaya 0 0 0 0 0 \
        -eI s_Nd 0 0 N_altai 0 0 0 0 \
        -eI s_D 0 0 0 N_denisova 0 0 0 \
        -ej t_Nd 3 1 \
        -ej t_Nc 2 1 \
        -ej t_ND 4 1 \
        -ej t_XY 6 5 -en 0.9482758620 5 4.8855 \
        -eps tm_NH 5 1 m_NH  \
        -ej t_HND 5 1 \
        -ej 86.2069 7 1 "
 	size_changes = list(n1="-n 1 0.850727 -en 0.468313 1 0.382459 -en 0.492840 1 0.968627 -en 0.522495 1 1.444230 -en 0.558351 1 1.821777 -en 0.601704 1 2.150016 -en 0.654122 1 2.354599 -en 0.717501 1 2.402091 -en 0.794132 1 2.307275 -en 0.886785 1 2.233972 -en 0.998812 1 2.226688 -en 1.134263 1 2.173059 -en 1.298035 1 2.073468 -en 1.496052 1 2.046063 -en 1.735472 1 2.151168 -en 2.024954 1 2.398446 -en 2.374965 1 2.795166 -en 2.798161 1 3.491499 -en 3.309844 1 4.828223 -en 3.928516 1 7.112816 -en 4.676550 1 10.071596 -en 5.580992 1 12.661877 -en 6.674547 1 14.233531 -en 7.996757 1 15.002602 -en 9.595433 1 14.984429 -en 11.528383 1 13.741670 -en 13.865499 1 10.257816 -en 20.107938 1 38.177939",
	n2="-n 2 0.175929 -en 0.719896 2 0.358807 -en 0.740175 2 0.745115 -en 0.764831 2 1.114071 -en 0.794809 2 1.550811 -en 0.831258 2 2.052445 -en 0.875575 2 2.445103",
	n3="-n 3 0.371086 -en 1.157749 3 0.510736 ",
	n4="-n 4 1.848249 -en 0.730929 4 1.235125 -en 0.758354 4 2.334873 -en 0.791521 4 2.836937 -en 0.831632 4 2.296415 -en 0.880141 4 1.598539 -en 0.938806 4 1.349040 -en 1.009754 4 1.468776 -en 1.095555 4 1.786108 -en 1.199321 4 2.103044 -en 1.324811 4 2.359405 -en 1.476575 4 2.670390 -en 1.660112 4 3.010718 -en 1.882077 4 3.203307 -en 2.150513 4 3.175625 -en 2.475150 4 3.033216 -en 2.867755 4 2.961124 -en 3.342559 4 3.189923 ",
	n5="-n 5 130 -en 0.00711 5 129.09 -en 0.0074 5 129.01 -en 0.00771 5 128.77 -en 0.00803 5 128.35 -en 0.00836 5 127.75 -en 0.0087 5 126.98 -en 0.00906 5 126.02 -en 0.00943 5 124.87 -en 0.00982 5 123.53 -en 0.01022 5 122.01 -en 0.01064 5 120.3 -en 0.01108 5 118.4 -en 0.01153 5 116.32 -en 0.01201 5 114.07 -en 0.0125 5 111.64 -en 0.01301 5 109.05 -en 0.01355 5 106.3 -en 0.01411 5 103.4 -en 0.01469 5 100.39 -en 0.01529 5 97.29 -en 0.01592 5 94.15 -en 0.01657 5 91 -en 0.01725 5 87.88 -en 0.01796 5 84.8 -en 0.0187 5 81.81 -en 0.01947 5 78.9 -en 0.02027 5 76.12 -en 0.0211 5 73.46 -en 0.02197 5 70.94 -en 0.02287 5 68.57 -en 0.02381 5 66.36 -en 0.02479 5 64.32 -en 0.02581 5 62.45 -en 0.02687 5 60.76 -en 0.02798 5 59.25 -en 0.02913 5 57.89 -en 0.03033 5 56.66 -en 0.03157 5 55.54 -en 0.03287 5 54.49 -en 0.03422 5 53.5 -en 0.03563 5 52.55 -en 0.03709 5 51.62 -en 0.03862 5 50.7 -en 0.0402 5 49.75 -en 0.04186 5 48.78 -en 0.04358 5 47.77 -en 0.04537 5 46.7 -en 0.04723 5 45.57 -en 0.04918 5 44.36 -en 0.0512 5 43.07 -en 0.0533 5 41.69 -en 0.05549 5 40.24 -en 0.05777 5 38.72 -en 0.06015 5 37.16 -en 0.06262 5 35.57 -en 0.0652 5 33.98 -en 0.06788 5 32.4 -en 0.07067 5 30.83 -en 0.07357 5 29.3 -en 0.07659 5 27.81 -en 0.07974 5 26.37 -en 0.08302 5 24.99 -en 0.08643 5 23.66 -en 0.08999 5 22.4 -en 0.09368 5 21.21 -en 0.09754 5 20.08 -en 0.10155 5 19.02 -en 0.10572 5 18.03 -en 0.11007 5 17.1 -en 0.11459 5 16.23 -en 0.1193 5 15.4 -en 0.1242 5 14.62 -en 0.12931 5 13.88 -en 0.13462 5 13.17 -en 0.14016 5 12.49 -en 0.14592 5 11.84 -en 0.15192 5 11.21 -en 0.15816 5 10.61 -en 0.16466 5 10.02 -en 0.17143 5 9.46 -en 0.17848 5 8.91 -en 0.18582 5 8.38 -en 0.19345 5 7.87 -en 0.20141 5 7.37 -en 0.20969 5 6.89 -en 0.2183 5 6.43 -en 0.22728 5 6 -en 0.23662 5 5.59 -en 0.24635 5 5.21 -en 0.25647 5 4.86 -en 0.26702 5 4.54 -en 0.27799 5 4.24 -en 0.28942 5 3.98 -en 0.30132 5 3.74 -en 0.3137 5 3.53 -en 0.3266 5 3.34 -en 0.34002 5 3.18 -en 0.354 5 3.04 -en 0.36855 5 2.92 -en 0.3837 5 2.83 -en 0.39947 5 2.75 -en 0.41589 5 2.7 -en 0.43299 5 2.67 -en 0.45079 5 2.65 -en 0.46932 5 2.65 -en 0.48861 5 2.67 -en 0.50869 5 2.7 -en 0.5296 5 2.75 -en 0.55137 5 2.81 -en 0.57404 5 2.89 -en 0.59763 5 2.99 -en 0.6222 5 3.1 -en 0.64778 5 3.24 -en 0.6744 5 3.39 -en 0.70213 5 3.56 -en 0.73099 5 3.75 -en 0.76104 5 3.97 -en 0.79232 5 4.21 -en 0.82489 5 4.48 -en 0.8588 5 4.78 -en 0.8941 5 5.11 -en 0.93085 5 5.47 -en 0.96911 5 5.86 -en 1.00895 5 6.28 -en 1.05042 5 6.74 -en 1.0936 5 7.23 -en 1.13855 5 7.75 -en 1.18536 5 8.31 -en 1.23408 5 8.89 -en 1.28481 5 9.5 -en 1.33762 5 10.13 -en 1.39261 5 10.78 -en 1.44985 5 11.44 -en 1.50945 5 12.11 -en 1.5715 5 12.76 -en 1.63609 5 13.41 -en 1.70335 5 14.05 -en 1.77336 5 14.66 -en 1.84626 5 15.24 -en 1.92215 5 15.79 -en 2.00116 5 16.29 -en 2.08342 5 16.75 -en 2.16906 5 17.16 -en 2.25823 5 17.5 -en 2.35105 5 17.79 -en 2.44769 5 18 -en 2.54831 5 18.15 -en 2.65306 5 18.22 -en 2.76212 5 18.21 -en 2.87565 5 18.13 -en 2.99386 5 17.97 -en 3.11693 5 17.74 -en 3.24505 5 17.45 -en 3.37844 5 17.12 -en 3.51732 5 16.75 -en 3.6619 5 16.36 -en 3.81242 5 15.95 -en 3.96914 5 15.55 -en 4.13229 5 15.14 -en 4.30215 5 14.75 -en 4.479 5 14.38 -en 4.66311 5 14.03 -en 4.85479 5 13.71 -en 5.05435 5 13.43 -en 5.26212 5 13.19 -en 5.47842 5 13 -en 5.70361 5 12.85 -en 5.93807 5 12.76 -en 6.18216 5 12.71 -en 6.43628 5 12.7 -en 6.70085 5 12.72 -en 6.97629 5 12.76 -en 7.26306 5 12.82 -en 7.56161 5 12.9 -en 7.87244 5 12.97 -en 8.19604 5 13.04 -en 8.53295 5 13.1 -en 8.8837 5 13.14 -en 9.24888 5 13.16 -en 9.62906 5 13.14 -en 10.02487 5 13.09 -en 10.43695 5 12.99 -en 10.86597 5 12.84 -en 11.31263 5 12.63 -en 11.77764 5 12.38 -en 12.26178 5 12.08 -en 12.76581 5 11.75 -en 13.29056 5 11.41 -en 13.83688 5 11.06 -en 14.40565 5 10.7 -en 14.99781 5 10.36 -en 15.61431 5 10.04 -en 16.25615 5 9.73 -en 16.92437 5 9.45 -en 17.62006 5 9.21 -en 18.34435 5 8.99 -en 19.09841 5 8.82 -en 19.88347 5 8.69 -en 20.7008 5 8.61 -en 21.55172 5 8.59",
	n6="-n 6 310 -en 0.00711 6 307.75 -en 0.0074 6 306.86 -en 0.00771 6 304.29 -en 0.00803 6 300.15 -en 0.00836 6 294.62 -en 0.0087 6 287.85 -en 0.00906 6 280.03 -en 0.00943 6 271.34 -en 0.00982 6 261.96 -en 0.01022 6 252.04 -en 0.01064 6 241.77 -en 0.01108 6 231.28 -en 0.01153 6 220.7 -en 0.01201 6 210.17 -en 0.0125 6 199.77 -en 0.01301 6 189.61 -en 0.01355 6 179.75 -en 0.01411 6 170.26 -en 0.01469 6 161.17 -en 0.01529 6 152.5 -en 0.01592 6 144.25 -en 0.01657 6 136.43 -en 0.01725 6 129.04 -en 0.01796 6 122.07 -en 0.0187 6 115.53 -en 0.01947 6 109.39 -en 0.02027 6 103.65 -en 0.0211 6 98.3 -en 0.02197 6 93.32 -en 0.02287 6 88.7 -en 0.02381 6 84.42 -en 0.02479 6 80.47 -en 0.02581 6 76.84 -en 0.02687 6 73.5 -en 0.02798 6 70.45 -en 0.02913 6 67.64 -en 0.03033 6 65.05 -en 0.03157 6 62.64 -en 0.03287 6 60.37 -en 0.03422 6 58.22 -en 0.03563 6 56.18 -en 0.03709 6 54.21 -en 0.03862 6 52.3 -en 0.0402 6 50.44 -en 0.04186 6 48.61 -en 0.04358 6 46.8 -en 0.04537 6 45 -en 0.04723 6 43.2 -en 0.04918 6 41.4 -en 0.0512 6 39.58 -en 0.0533 6 37.76 -en 0.05549 6 35.93 -en 0.05777 6 34.12 -en 0.06015 6 32.35 -en 0.06262 6 30.63 -en 0.0652 6 28.97 -en 0.06788 6 27.39 -en 0.07067 6 25.9 -en 0.07357 6 24.49 -en 0.07659 6 23.18 -en 0.07974 6 21.97 -en 0.08302 6 20.85 -en 0.08643 6 19.82 -en 0.08999 6 18.89 -en 0.09368 6 18.06 -en 0.09754 6 17.31 -en 0.10155 6 16.66 -en 0.10572 6 16.1 -en 0.11007 6 15.61 -en 0.11459 6 15.2 -en 0.1193 6 14.84 -en 0.1242 6 14.54 -en 0.12931 6 14.28 -en 0.13462 6 14.06 -en 0.14016 6 13.87 -en 0.14592 6 13.7 -en 0.15192 6 13.56 -en 0.15816 6 13.44 -en 0.16466 6 13.33 -en 0.17143 6 13.22 -en 0.17848 6 13.12 -en 0.18582 6 13.02 -en 0.19345 6 12.91 -en 0.20141 6 12.8 -en 0.20969 6 12.67 -en 0.2183 6 12.54 -en 0.22728 6 12.4 -en 0.23662 6 12.25 -en 0.24635 6 12.11 -en 0.25647 6 11.96 -en 0.26702 6 11.82 -en 0.27799 6 11.68 -en 0.28942 6 11.55 -en 0.30132 6 11.43 -en 0.3137 6 11.32 -en 0.3266 6 11.22 -en 0.34002 6 11.14 -en 0.354 6 11.07 -en 0.36855 6 11.02 -en 0.3837 6 10.99 -en 0.39947 6 10.98 -en 0.41589 6 10.99 -en 0.43299 6 11.02 -en 0.45079 6 11.07 -en 0.46932 6 11.14 -en 0.48861 6 11.24 -en 0.50869 6 11.35 -en 0.5296 6 11.47 -en 0.55137 6 11.62 -en 0.57404 6 11.78 -en 0.59763 6 11.96 -en 0.6222 6 12.16 -en 0.64778 6 12.37 -en 0.6744 6 12.6 -en 0.70213 6 12.84 -en 0.73099 6 13.1 -en 0.76104 6 13.37 -en 0.79232 6 13.66 -en 0.82489 6 13.96 -en 0.8588 6 14.28 -en 0.8941 6 14.6 -en 0.93085 6 14.94")
	main1 = sub("N_denisova", N_denisova, sub("N_altai", N_altai, sub("N_chagyrskaya", N_chagyrskaya, sub("N_vindija", N_vindija, main1))))
    for (i in names(par_to_sim)) { main1=gsub(paste(i,""), paste(par_to_sim[[i]][1],""), as.character(main1)) }
    cmd = paste(main, main1, test_pop, paste(unlist(size_changes),collapse=" "), cmd_append)
    if (run) {
        library(scrm)
		return(scrm(cmd)[[1]])
	} else {
		return(cmd)
	}
}



## A more complex model (but Ne changing only at key events) including two more events of admixture:
## 1) Superarchiac into Denisova
## 2) Human into Neandertals
## Notice that gene flow from Neandertals to Humans is only into Africans, which can be neglected.
sim_scrm_afr = function(nsims , t_pop=NULL, split_pop='v', basepairs=1e6, RHO=52, THETA=58, BASEPAIRS=1e6, rescale=TRUE, n_test=1, n_vin=1, n_cha=1, n_alt=1, n_den=1, n_afr=1, n_chimp=0, params=NULL, run=TRUE, cmd_append=NULL) {
	if (is.null(t_pop)) {
		stop("t_pop should be specified.")
	}
	# Fixed parameters from recent momi estimates using 5 Yorubas.
	# Also t_x = 300kya, a whatever fixed paparmeters
	par_to_sim = list(t_HND=4.4138, t_ND=3.7759, t_Nd=1.0862, t_Nc=0.6897,  t_sA=21.4574, t_iH=3.7672, n_ND=3, n_Nvcd=1.964, n_Nd=0.536, n_Nc=0.303, n_Nv=1.107, n_D=2.068, n_HND=15.428, n_Nvc=2.515, n_Hy=23.421, n_sA=15, n_iH=10, s_Nv=0.4741, s_Nc=0.6293, s_Nd=0.9483, s_D=0.5776, tm_SD=2.5534, tm_NH=0.2155, tm_HN=1.6262, m_NH=0.9979, m_HN=.9273, m_SD=0.9794)
	if (!is.null(params)) {
		if( sum(names(params) %in% names(par_to_sim)) == length(params) ) {
			par_to_sim[names(params)] = params
		} else {
			stop("all names of params should match.")
		}
	}
	require("scrm")
	pops = c("v", "c", "a", "d", "x", "i", "s", "o")
	# The order of the population are vindija (1), chagyrskaya (2), altai (3), denisova (4), africa (5), introgressing human (6), ghost|super-archaic (7), chimp (8).
	# Notice that this order is different from the order of the samples in the output of the simulation, where unkown human, super-archaic and chimp will not be present.
	test_pop = rep(0, length(pops))
	test_pop[which(split_pop == pops)] = n_test
	test_pop = paste("-eI", t_pop, paste(test_pop, collapse=" "))
	## Rescale the values RHO and THETA depending on the ratio basepairs/BASEPAIRS
	if (rescale) {
		RHO = RHO*basepairs/BASEPAIRS
		THETA = THETA*basepairs/BASEPAIRS
	}
	n_tot = sum(as.numeric(c(n_test, n_vin, n_cha, n_alt, n_den, n_afr, n_chimp)))
	main = sub("n_chimp", n_chimp, sub("n_afr", n_afr, paste(n_tot, nsims, "-t", THETA, "-r", RHO, basepairs, "-I 8 0 0 0 0 n_afr 0 0 n_chimp -l 100r -SC abs ")))
	main1 = "-eI s_Nv n_vin 0 0 0 0 0 0 0 \
             -eI s_Nc 0 n_cha 0 0 0 0 0 0 \
             -eI s_Nd 0 0 n_alt 0 0 0 0 0 \
             -eI s_D 0 0 0 n_den 0 0 0 0 \
             -ej t_Nd 3 1 \
             -ej t_Nc 2 1 \
             -ej t_ND 4 1 \
             -ej t_iH 6 5 \
             -ej t_HND 5 1 -en t_HND 1 n_HND \
             -ej t_sA 7 1 \
             -ej 86.207 8 1 \
             -eps tm_NH 5 1 m_NH \
             -eps tm_SD 4 7 m_SD \
             -eps tm_HN 1 6 m_HN \
             -n 1 n_Nv -en t_Nc 1 n_Nvc -en t_Nd 1 n_Nvcd -en t_ND 1 n_ND  \
             -n 2 n_Nc -n 3 n_Nd -n 4 n_D -n 5 n_Hy -n 6 n_iH -n 7 n_sA "
	main1 = sub("n_den", n_den, sub("n_alt", n_alt, sub("n_cha", n_cha, sub("n_vin", n_vin, main1))))
	for (i in names(par_to_sim)) { main1=gsub(paste(i,""), paste(par_to_sim[[i]][1],""), as.character(main1)) }
	cmd = paste(main, main1, cmd_append, test_pop)
	if (run) {
		return(scrm(cmd)[[1]])
	} else {
		return(cmd)
	}
}

## Including one or two more events of Neandertals to Denisova admixture:
## tm_ND from the Altai branch
## tm_ND1 from the Chagyrskaya branch
sim_scrm_afr_mND = function(nsims , t_pop=NULL, split_pop='v', basepairs=1e6, RHO=52, THETA=58, BASEPAIRS=1e6, rescale=TRUE, n_test=1, n_vin=1, n_cha=1, n_alt=1, n_den=1, n_afr=1, params=NULL, run=TRUE) {
	if (is.null(t_pop)) {
		stop("t_pop should be specified.")
	}
	# Fixed parameters from recent momi estimates using 5 Yorubas.
	# Also t_x = 300kya, a whatever fixed paparmeters
	par_to_sim = list(t_HND=4.4138, t_ND=3.7759, t_Nd=1.0922, t_Nc=0.6897,  t_sA=21.4574, t_iH=3.7672, n_ND=3, n_Nvcd=1.964, n_Nd=0.536, n_Nc=0.303, n_Nv=1.107, n_D=2.068, n_HND=15.428, n_Nvc=2.515, n_Hy=23.421, n_sA=15, n_iH=10, s_Nv=0.4741, s_Nc=0.6293, s_Nd=0.9483, s_D=0.5776, tm_SD=2.5534, tm_NH=0.2155, tm_HN=1.6262, m_NH=0.9979, m_HN=.9273, m_SD=0.9794, tm_ND=1.0991, m_ND=0.9794, tm_ND1=NULL,  m_ND1=NULL)
	if (!is.null(params)) {
		if( sum(names(params) %in% names(par_to_sim)) == length(params) ) {
			par_to_sim[names(params)] = params
		} else {
			stop("all names of params should match.")
		}
	}
	require("scrm")
	pops = c("v", "c", "a", "d", "x", "i", "s", "o")
	# The order of the population are vindija (1), chagyrskaya (2), altai (3), denisova (4), africa (5), introgressing human (6), ghost|super-archaic (7), chimp (8).
	# Notice that this order is different from the order of the samples in the output of the simulation, where unkown human, super-archaic and chimp will not be present.
	test_pop = rep(0, length(pops))
	test_pop[which(split_pop == pops)] = n_test
	test_pop = paste("-eI", t_pop, paste(test_pop, collapse=" "))
	## Rescale the values RHO and THETA depending on the ratio basepairs/BASEPAIRS
	if (rescale) {
		RHO = RHO*basepairs/BASEPAIRS
		THETA = THETA*basepairs/BASEPAIRS
	}
	n_tot = sum(as.numeric(c(n_test, n_vin, n_cha, n_alt, n_den, n_afr)))
	main = sub("n_afr", n_afr, paste(n_tot, nsims, "-t", THETA, "-r", RHO, basepairs, "-I 8 0 0 0 0 n_afr 0 0 0 -l 100r -SC abs "))
	main1 = "-eI s_Nv n_vin 0 0 0 0 0 0 0 \
             -eI s_Nc 0 n_cha 0 0 0 0 0 0 \
             -eI s_Nd 0 0 n_alt 0 0 0 0 0 \
             -eI s_D 0 0 0 n_den 0 0 0 0 \
             -ej t_Nd 3 1 \
             -ej t_Nc 2 1 \
             -ej t_ND 4 1 \
             -ej t_iH 6 5 \
             -ej t_HND 5 1 -en t_HND 1 n_HND \
             -ej t_sA 7 1 \
             -ej 86.207 8 1 \
             -eps tm_NH 5 1 m_NH \
             -eps tm_SD 4 7 m_SD \
             -eps tm_HN 1 6 m_HN \
             -n 1 n_Nv -en t_Nc 1 n_Nvc -en t_Nd 1 n_Nvcd -en t_ND 1 n_ND  \
             -n 2 n_Nc -n 3 n_Nd -n 4 n_D -n 5 n_Hy -n 6 n_iH -n 7 n_sA "
	main1 = sub("n_den", n_den, sub("n_alt", n_alt, sub("n_cha", n_cha, sub("n_vin", n_vin, main1))))
	for (i in names(par_to_sim)) { main1=gsub(paste(i,""), paste(par_to_sim[[i]][1],""), as.character(main1)) }
	migr_nd = ""
	if(is.numeric(par_to_sim$m_ND) & is.numeric(par_to_sim$tm_ND) & par_to_sim$tm_ND < par_to_sim$t_Nd ) {
		migr_nd = paste("-eps", par_to_sim$tm_ND, "4 3", par_to_sim$m_ND)
	}
	if(sum(is.numeric(par_to_sim$m_ND1)) + sum(par_to_sim$tm_ND1) + sum(par_to_sim$tm_ND1 < par_to_sim$t_Nc ) == 3 ) {
		 migr_nd = paste(migr_nd, paste("-eps", par_to_sim$tm_ND1, "4 2", par_to_sim$m_ND1))
	}
	cmd = paste(main, main1, migr_nd, test_pop )
	if (run) {
		return(scrm(cmd)[[1]])
	} else {
		return(cmd)
	}
}


## Two African populations
sim_scrm_afr2 = function(nsims , t_pop=NULL, split_pop='v', basepairs=1e6, RHO=52, THETA=58, BASEPAIRS=1e6, rescale=TRUE, n_test=1, n_vin=1, n_cha=1, n_alt=1, n_den=1, n_afr=2, n_mbu=2, params=NULL, run=TRUE) {
	if (is.null(t_pop)) {
		stop("t_pop should be specified.")
	}
	# Fixed parameters from recent momi estimates using 5 Yorubas.
	# Also t_x = 300kya, a whatever fixed paparmeters
	par_to_sim = list(t_HND=4.4138, t_ND=3.7759, t_Nd=1.0862, t_Nc=0.6897,  t_sA=21.4574, t_iH=3.7672, t_Hm=1.7298, n_ND=3, n_Nvcd=1.964, n_Nd=0.536, n_Nc=0.303, n_Nv=1.107, n_D=2.068, n_HND=15.428, n_Nvc=2.515, n_Hy=23.421, n_Hm=32.6, n_Hym=1.231, n_sA=15, n_iH=10, s_Nv=0.4741, s_Nc=0.6293, s_Nd=0.9483, s_D=0.5776, tm_SD=2.5534, tm_ND=1.2082, tm_NH=0.2155, tm_HN=1.6262, m_NH=0.9979, m_ND=0.9794, m_HN=.9273, m_SD=0.9794)
	if (!is.null(params)) {
		if( sum(names(params) %in% names(par_to_sim)) == length(params) ) {
			par_to_sim[names(params)] = params
		} else {
			stop("all names of params should match.")
		}
	}
	require("scrm")
	pops = c("v", "c", "a", "d", "y", "x", "i", "s", "o")
	# The order of the population are vindija (1), chagyrskaya (2), altai (3), denisova (4), africa (5), africa (6), introgressing human (7), ghost|super-archaic (8), chimp (9).
	# Notice that this order is different from the order of the samples in the output of the simulation, where unkown human, super-archaic and chimp will not be present.
	test_pop = rep(0, length(pops))
	test_pop[which(split_pop == pops)] = n_test
	test_pop = paste("-eI", t_pop, paste(test_pop, collapse=" "))
	## Rescale the values RHO and THETA depending on the ratio basepairs/BASEPAIRS
	if (rescale) {
		RHO = RHO*basepairs/BASEPAIRS
		THETA = THETA*basepairs/BASEPAIRS
	}
	n_tot = sum(as.numeric(c(n_test, n_vin, n_cha, n_alt, n_den, n_afr, n_mbu)))
	main = sub("n_mbu",n_mbu, sub("n_afr", n_afr, paste(n_tot, nsims, "-t", THETA, "-r", RHO, basepairs, "-I 9 0 0 0 0 n_afr n_mbu 0 0 0 -l 100r -SC abs ")))
	main1 = "-eI s_Nv n_vin 0 0 0 0 0 0 0 0 \
             -eI s_Nc 0 n_cha 0 0 0 0 0 0 0 \
             -eI s_Nd 0 0 n_alt 0 0 0 0 0 0 \
             -eI s_D  0 0 0 n_den 0 0 0 0 0 \
             -ej t_Nd 3 1 \
             -ej t_Nc 2 1 \
             -ej t_ND 4 1 \
             -ej t_Hm 6 5 -en t_Hm 1 n_Hym \
             -ej t_iH 7 5 \
             -ej t_HND 5 1 -en t_HND 1 n_HND \
             -ej t_sA 8 1 \
             -ej 86.207 9 1 \
             -eps tm_NH 5 1 m_NH \
             -eps tm_ND 4 1 m_ND \
             -eps tm_SD 4 7 m_SD \
             -eps tm_HN 1 6 m_HN \
             -n 1 n_Nv -en t_Nc 1 n_Nvc -en t_Nd 1 n_Nvcd -en t_ND 1 n_ND  \
             -n 2 n_Nc -n 3 n_Nd -n 4 n_D -n 8 n_sA \
             -n 5 n_Hy -n 6 n_Hm -n 7 n_iH "
	main1 = sub("n_den", n_den, sub("n_alt", n_alt, sub("n_cha", n_cha, sub("n_vin", n_vin, main1))))
	for (i in names(par_to_sim)) { main1=gsub(paste(i,""), paste(par_to_sim[[i]][1],""), as.character(main1)) }
	cmd = paste(main, main1, test_pop)
	if (run) {
		return(scrm(cmd)[[1]])
	} else {
		return(cmd)
	}
}


## Function to run simulations with more complex model, although Ne does not change piecewise.
## These are the admixture events, and each of them can happen more than once. Set the timing to NULL or the intesity to 1 to make them not happaning.
## 1) Superarchiac into Denisova
## 2) Human into Neandertals
## 3) Neadertals into Denisovans.
##   So far only from Altai or all Nea if the event is older than the split time of Altai, t_Nd.
## 4) Neadertals into Humans

# t_pop=NULL; split_pop='v'; basepairs=1e6; RHO=52; THETA=58; BASEPAIRS=1e6; rescale=TRUE; n_test=1; n_vin=1; n_cha=1; n_alt=1; n_den=1; n_afr=1; params=NULL; run=TRUE

## Notice that gene flow from Neandertals to Humans is only into Africans, which can be neglected.
sim_scrm = function(nsims , t_pop=NULL, split_pop='v', basepairs=1e6, RHO=52, THETA=58, BASEPAIRS=1e6, rescale=TRUE, n_test=1, n_vin=1, n_cha=1, n_alt=1, n_den=1, n_afr=1, n_chimp=0, params=NULL, run=TRUE, cmd_append=NULL) {
	# Set the parameters to NULL to cancel them.
	if (is.null(t_pop)) {
		stop("t_pop should be specified.")
	}
	par_codes = list(sizes="-n 1 n_Nv -en t_Nc 1 n_Nvc -en t_Nd 1 n_Nvcd -en t_ND 1 n_ND -n 2 n_Nc -n 3 n_Nd -n 4 n_D -n 5 n_Hy -n 6 n_iH -n 7 n_sA",
	    s_Nv="-eI s_Nv n_vin 0 0 0 0 0 0 0",
		s_Nc="-eI s_Nc 0 n_cha 0 0 0 0 0 0",
		s_Nd="-eI s_Nd 0 0 n_alt 0 0 0 0 0",
		s_D="-eI s_D 0 0 0 n_den 0 0 0 0",
		t_Nd="-ej t_Nd 3 1",
		t_Nc="-ej t_Nc 2 1",
		t_ND="-ej t_ND 4 1",
		t_iH="-ej t_iH 6 5",
		t_HND="-ej t_HND 5 1 -en t_HND 1 n_HND",
		t_sA="-ej t_sA 7 1",
		t_chimp="-ej t_chimp 8 1",
		tm_NH="-eps tm_NH 5 1 m_NH",
		tm_SD="-eps tm_SD 4 7 m_SD",
		tm_HN="-eps tm_HN 1 6 m_HN",
		tm_ND="-eps tm_ND 4 3 m_ND"
	)
	# Fixed parameters from recent momi estimates using 5 Yorubas.
	# Also t_x = 300kya, a whatever fixed paparmeters
	par_to_sim = list(t_HND=4.4138, t_ND=3.7759, t_Nd=1.0922, t_Nc=0.6897,  t_sA=21.4574, t_iH=3.7672, n_ND=3, n_Nvcd=1.964, n_Nd=0.536, n_Nc=0.303, n_Nv=1.107, n_D=2.068, n_HND=15.428, n_Nvc=2.515, n_Hy=23.421, n_sA=15, n_iH=10, s_Nv=0.4741, s_Nc=0.6293, s_Nd=0.9483, s_D=0.5776, tm_SD=2.5534, tm_NH=0.2155, tm_HN=1.6262, tm_ND=1.0991, m_NH=0.9979, m_HN=0.9273, m_SD=0.9794, m_ND=0.9794, t_chimp=86.207)
	if (!is.null(params)) {
		if( sum(names(params) %in% names(par_to_sim)) == length(params) ) {
			par_to_sim[names(params)] = params
		} else {
			stop("all names of params should match.")
		}
	}
	pops = c("v", "c", "a", "d", "x", "i", "s", "o")
	# The order of the population are vindija (1), chagyrskaya (2), altai (3), denisova (4), africa (5), introgressing human (6), ghost|super-archaic (7), chimp (8).
	# Notice that this order is different from the order of the samples in the output of the simulation, where unkown human, super-archaic and chimp will not be present.
	test_pop = rep(0, length(pops))
	test_pop[which(split_pop == pops)] = n_test
	test_pop = paste("-eI", t_pop, paste(test_pop, collapse=" "))
	## Rescale the values RHO and THETA depending on the ratio basepairs/BASEPAIRS
	if (rescale) {
		RHO = RHO*basepairs/BASEPAIRS
		THETA = THETA*basepairs/BASEPAIRS
	}
	n_tot = sum(as.numeric(c(n_test, n_vin, n_cha, n_alt, n_den, n_afr, n_chimp)))
	# make sep="\n" if you want to debug and look at the code in nicer printing layout in newlines
	sep="\n"
	main = paste(sub("n_chimp", n_chimp, sub("n_afr", n_afr, paste(n_tot, nsims, "-t", THETA, "-r", RHO, basepairs, "-I 8 0 0 0 0 n_afr 0 0 n_chimp -l 100r -SC abs "))), sep, par_codes$sizes, sep)
	for (i in names(par_codes)[-1]) {
		cc = par_codes[[i]]
		for (j in 1:length(par_to_sim[[i]])) {
			p = par_to_sim[[i]][j]
			if(!is.null(p)) {
				# Check consistency with the introgressing Nea into Den
				if (i == "tm_ND" & p >  par_to_sim[["t_Nd"]]) {
					cc="-eps tm_ND 4 1 m_ND"
				}
				if (i == "tm_SD" & p >  par_to_sim[["t_ND"]]) {
					cc="-eps tm_SD 1 7 m_SD"
				}
				main=paste(main, gsub(paste(i, ""), paste(p, ""), cc), "\n")
			}
		}
	}
	par_names=names(par_to_sim)
	for (m in grep("^m_", par_names )) {
		i=par_names[m]
		for(j in 1:length(par_to_sim[[i]])) {
			main=sub(i, par_to_sim[[i]][j], main)
		}
	# print(par_to_sim[[par_names[m]]])
	}
	for(n in grep("^n_", par_names )) {
		i=par_names[n]
		main=gsub(paste(i,""), paste(par_to_sim[[i]][1], ""), main)
	}
	for(i in grep("^t_", par_names )) {
		i=par_names[i]
		main=gsub(paste(i,""), paste(par_to_sim[[i]][1], ""), main)
	}
	main = sub("n_den", n_den, sub("n_alt", n_alt, sub("n_cha", n_cha, sub("n_vin", n_vin, main))))
	cmd = paste(main, test_pop)
	if (run) {
		require(scrm)
		return(scrm(cmd)[[1]])
	} else {
		return(cmd)
	}
}


##########################################################################################
## Generate a table like the lineage assignment or class of sites from scrm simulations ##

mktab <- function(sim_data, human="africa", pop_ids=list(europe=1:2, africa=3:4, vindija=5:6, chagyrskaya=7:8, altai=9:10, denisova=11:12, test=13:14), contamination = 0, pop_cont_source="europe", pop_cont_sink="test", all_combinations = NULL) {
	if(is.null(all_combinations)) {
	all_combinations <- c("altai", "altai-chagyrskaya", "altai-chagyrskaya-denisova", "altai-denisova", "altai-vindija", "altai-vindija-denisova", "chagyrskaya", "chagyrskaya-denisova", "denisova", "human", "human-altai", "human-altai-chagyrskaya", "human-altai-chagyrskaya-denisova", "human-altai-denisova", "human-altai-vindija", "human-altai-vindija-denisova", "human-chagyrskaya", "human-chagyrskaya-denisova", "human-denisova", "human-neandertal", "human-vindija", "human-vindija-chagyrskaya", "human-vindija-chagyrskaya-denisova", "human-vindija-denisova", "neandertal", "neandertal-denisova", "vindija", "vindija-chagyrskaya", "vindija-chagyrskaya-denisova", "vindija-denisova")
	}
    popN <- sapply(pop_ids, length)
    n <- sum(popN)
	## Traspose the data since using scrm within R was different than using scrm with the shell
	## This is done because it's faster to have many rows rather than many columns.
	sim_data <- lapply(sim_data, t)
    ## Substitute "human" with africa or europe depending for which population was used for the lineage assignments
    all_combinations <- sub("human", human, sub("neandertal", "altai-vindija-chagyrskaya", all_combinations))
    ## Concatenate the simulations in order to match the real data
    ## Add contamination by replacing the sink population with the source one in a given proportion of simulations.
    ## Notice that this is the cheap way to do it because I mix simulations and not portion of sequences/individuals.
    if(contamination > 0) {
        i <- sort(sample(1:length(sim_data), round(length(sim_data)*as.numeric(contamination))))
        for (ii in i ) {
            x=sim_data[[ii]];
            x[, pop_ids[[pop_cont_sink]]] <- x[,pop_ids[[pop_cont_source]]]
            sim_data[[ii]] <- x
        }
    }
    sim_data = do.call("rbind",sim_data)
    pids <- c(sapply(pop_ids[names(pop_ids) %in% all_combinations], function(i) sample(i,1)), test=sample(pop_ids[["test"]],1))
    ## Monomorphic sites to be removed (it can be speeded up)
    mono1 <- which(apply(sim_data[,pids], 1, function(m) sum(m == "0") ) == length(pids)) # monomorphic ancestral
    mono2 <- which(apply(sim_data[,pids], 1, function(m) sum(m == "1") ) == length(pids)) # monomorphic derived
    x <- matrix(as.numeric(sim_data[-c(mono1, mono2),pids]),ncol=length(pids)); colnames(x) <- names(pids)
    ## one european and african chromosome
    if(human == "africa") {
        h = pop_ids[["europe"]]
    } else {
        h <- pop_ids[["africa"]];
    }
    ids <- c(sample(pop_ids[[human]][!pop_ids[[human]] %in% pids[human] ],1), sample(h, 1) )
    x1 <- sim_data[-c(mono1, mono2), ids]
    rm(sim_data); for (i in 1:10) {g <- gc()}
    ## The informative sites ids given by "all_combinations" (an ugly code but efficient and fast)
	ids <- c()
    for (y in all_combinations) {
        i <- strsplit(y, split="-")[[1]]
        j <- colnames(x)[!colnames(x) %in% c(i,"test")]
        y1 <- rowSums(sim_data[, i, drop=FALSE])
        y2 <- rowSums(sim_data[, j, drop=FALSE])
        ids[[y]] <- which(y1 == length(i) & y2 == 0)
    }
    ## Get the count for one simulated haplotype of the test population.
    tab <- t(sapply(ids, function(i) if(length(i) > 0) { c(sum(x[i,"test"] == 0), sum(x[i,"test"] == 1))  } else { c(0,0) } ))
    ## Get the count for one random simulated european and african chromosome
    tabs <- do.call("cbind", lapply(1:2, function(s) t(sapply(ids, function(i) if(length(i) > 0) { c(sum(x1[i,s] == 0), sum(x1[i,s] == 1))  } else { c(0,0) } ))))
    return(cbind(tab, tabs))
}


## I should change these with one letters abbreviation like 'a', 'a-c', 'a-c-v-d', etc...
all_combinations <- ClassOfSites <- c("altai", "altai-chagyrskaya", "altai-chagyrskaya-denisova", "altai-denisova", "altai-vindija", "altai-vindija-denisova", "chagyrskaya", "chagyrskaya-denisova", "denisova", "human", "human-altai", "human-altai-chagyrskaya", "human-altai-chagyrskaya-denisova", "human-altai-denisova", "human-altai-vindija", "human-altai-vindija-denisova", "human-chagyrskaya", "human-chagyrskaya-denisova", "human-denisova", "human-neandertal", "human-vindija", "human-vindija-chagyrskaya", "human-vindija-chagyrskaya-denisova", "human-vindija-denisova", "neandertal", "neandertal-denisova", "vindija", "vindija-chagyrskaya", "vindija-chagyrskaya-denisova", "vindija-denisova")


###################################
###################################
## To read the input file(s).
read.lineages <- function(File) {
	# ClassOfSites must be present in the workspace (see above).
	a <- scan(File,what="",sep="\n",quiet=T)
	ids <- grep("#",a);	names(ids) <- sapply(strsplit(a[ids], split="\t"), function(x) sub("#","",x[1]))
	header = unlist(strsplit(gsub("f_", "", gsub("%", "", a[ids[1]])), split="\t"))
	if(length(a)-length(ids)-1 <= length(ClassOfSites)) {
		a = do.call("rbind", lapply(strsplit(a[-ids],split="\t"), function(x) x[1:3]));
		a = a[na.omit(match(ClassOfSites, a[,1])),]
		b = cbind(Ancestral=as.numeric(a[,2]), Derived=as.numeric(a[,3]))
		rownames(b) = a[,1]; a = list(b)
	} else {
		if (length(header) > 4) {
			cc = c(2:3,5:length(header)); names(cc) = c("Ancestral", "Derived", header[cc][-(1:2)])
		} else {
			cc = 2:3; names(cc) = c("Ancestral", "Derived")
		}
		a = lapply(1:length(ids), function(x) {if(x < length(ids)) {i1=(ids[x]+1); i2=(ids[x+1]-1) } else { i1=(ids[x]+1); i2=length(a) } ; b = do.call("rbind",strsplit(a[i1:i2], split="\t"));
		class_of_sites=b[,1]
		if(sum(ClassOfSites %in% class_of_sites) < length(ClassOfSites)) {
			b=rbind(b, matrix(c(ClassOfSites[!ClassOfSites %in% class_of_sites], rep(0, sum(!ClassOfSites %in% class_of_sites)*(ncol(b)-1) )),ncol=ncol(b)))
		}
		b = apply(b[na.omit(match(ClassOfSites, b[,1])),cc],2,as.numeric);
		dimnames(b) = list(ClassOfSites, names(cc)); b})
		names(a) = names(ids)
	}
    return(a)
}

# A slithly different version which takes into account the type of ancestry and different human sample references but it chooses only one between Mbuti (default) and Yoruba.
read_lineages <- function(File, anc=NULL, hum="m") {
	# ClassOfSites must be present in the workspace (see above).
	a <- scan(File, what="",sep="\n", quiet=T)
	ids <- grep("#",a);	names(ids) <- sapply(strsplit(a[ids], split="\t"), function(x) sub("#","",x[1]))
	header <- unlist(strsplit(gsub("f_", "", gsub("%", "", a[ids[1]])), split="\t"))[-1]
	x <- strsplit(a, split="\t")
	ClassOfSites <- unique(sort(sapply(x[-ids],function(y) y[1])))
	if(length(a)-length(ids)-1 <= length(ClassOfSites)) {
		# This should be modified
		z = do.call("rbind", strsplit(a[-ids],split="\t"));
		a = z[na.omit(match(ClassOfSites, z[,1])),]
		b = apply(a[,-1], 2, as.numeric)
		dimnames(b) = list(a[,1], header)
		a = list("#all_chrs"=b)
	} else {
		a = lapply(1:length(ids), function(j) {
			if(j < length(ids)) {
				i1=(ids[j]+1)
				i2=(ids[j+1]-1)
			} else {
				i1=(ids[j]+1)
				i2=length(a)
			} ; b = do.call("rbind",strsplit(a[i1:i2], split="\t")); class_of_sites=b[,1]
			if(sum(ClassOfSites %in% class_of_sites) < length(ClassOfSites)) {
				b = rbind(b, matrix(c(ClassOfSites[!ClassOfSites %in% class_of_sites], rep(0, sum(!ClassOfSites %in% class_of_sites)*(ncol(b)-1) )),ncol=ncol(b)));
				}
			b = apply(b[na.omit(match(ClassOfSites, b[,1])),-1][,c(T,T,F)], 2, as.numeric);
			dimnames(b) = list(ClassOfSites, header[c(T,T,F)]); b
		})
	}
	cc = do.call("rbind", strsplit(ClassOfSites,split="_"))
	if(ncol(cc) == 2) {
		colnames(cc) = c("Anc","Lin")
		if(!is.null(anc) & sum(cc[,1] %in% anc) > 1) {
			cc = subset(as.data.frame(cc), Anc %in% anc)
		}
		# Take only one human population (hum) as reference and remove the other (h) 
		h = c("m","y")[c("m","y") != hum] # the 
		cc = cbind(cc, Lin1=sub(hum, "x", sub(h, "", cc[,2])))
		# this does nothing with the older format.
		cc = cc[nchar(cc[,3]) > 0,]
		ClassOfSites = paste(cc[,1],cc[,2],sep="_")
	} else { # this is just to be compatible with the older format.
		cc = cbind(Anc=rep("", nrow(cc)), Lin=as.vector(cc), Lin1=as.vector(cc))
		ClassOfSites = cc[,"Lin"]
	}
	a = lapply(a, function(z) {y=z[ClassOfSites,]; z=t(sapply(unique(cc[,3]), function(i) colSums(z[ClassOfSites[cc[,3] == i],,drop=F]))); rownames(z) == unique(cc[,3]); z } )
	if(length(a) > 1 ) {
		names(a) = names(ids)
    }
	return(a)
}

# Function to simplify the ClassOfSites (or lineage assignment) table.
simplify_lineage = function(xdata, samples=c("human", "vindija", "denisova")) {
	all_comb = unlist(sapply(1:length(samples), function(i) apply(combn(samples,i),2, function(x) paste(x,collapse="-"))))
	all_comb_new = strsplit(all_comb, split="-")
	all_comb_old = strsplit(rownames(xdata), split="-")
	x1 = xdata
	out = matrix(0, ncol=ncol(xdata), nrow=length(all_comb))
	rownames(out) = all_comb
	for(i in length(all_comb):1) {
		ids = sapply(all_comb_old, function(x) sum(x %in% all_comb_new[[i]]))
		id = which(ids == length(all_comb_new[[i]]))
		out[i,] = colSums(x1[id,drop=FALSE,])
		x1 = x1[-id,]
		all_comb_old = all_comb_old[-id]
	}
	return(out)
}


## The function to calculate the likelihood considering human contamination
LogLikeCon <- function(obs, sim, con=0, weighted=TRUE, multinomial=TRUE) {
    ## sim:      must have 4 columns, 1st-2nd columns with Anc and Der counts in the test sample, 3rd-4th with counts in humans (Africa or Europe) to to be considered as contamiant.
    ## con:      the proportion of contamiantion to be considered for estimating the likelihood.
    ## weighted: is useful when multinomial=FALSE so that the binomial will be normalized by the proportion of sites.
    if(multinomial) {
        p <- c(sim[,1],sim[,2])/sum(sim[,1:2])
        ## Correct the counts for human contamination given by the 3rd-4th columns
        p1 <- c(sim[,3], sim[,4])/sum(sim[,3:4])
        p1a <- (1-con)*p+(con*p1);
        d <- log(dmultinom(c(obs[,1],obs[,2]), prob=p1a))
    } else {
        p <- sim[,2]/(sim[,1]+sim[,2])
        p1 <- sim[,4]/(sim[,3]+sim[,4])
        p1a <- (1-con)*p+(con*p1); p1a[p1a > 1 ] <- 1; p1a[p1a<0] <- 0
        if (weighted) {
            d <- sapply(1:nrow(obs), function(i) log(dbinom(obs[i,2], sum(obs[i,1:2]),prob=p1a[i])*sum(obs[i,1:2])))
        } else {
            d <- sapply(1:nrow(obs), function(i) log(dbinom(obs[i,2], sum(obs[i,1:2]),prob=p1a[i])))
        }
    }
    return(d)
}

## The function to calculate the likelihood integrating human contamination.
LikeCon <- function(obs, sim, con=0, Binomial=FALSE) {
    ## sim:       must have 4 columns, 1st-2nd columns with Anc and Der counts in the test sample, 3rd-4th with counts in humans (Africa or Europe) to be considered as source of contamination.
    ## con:       the proportion of contamiantion to be considered for estimating the likelihood.
    if (Binomial) {
		p <- sim[,2]/(sim[,1]+sim[,2])
        p_contaminant <- sim[,4]/(sim[,3]+sim[,4])
        p_corrected <- (1-con)*p+(con*p_contaminant)
        p_corrected[p_corrected > 1 ] <- 1
        p_corrected[p_corrected<0] <- 0
		d <- sum(sapply(1:nrow(obs), function(i) dbinom(obs[i,2], sum(obs[i,1:2]), prob=p_corrected[i], log=TRUE)))
    } else {
		p <- c(sim[,1],sim[,2])/sum(sim[,1:2])
		## Correct the counts for human contamination given by the 3rd-4th columns
		p_contaminant <- c(sim[,3], sim[,4])/sum(sim[,3:4])
		p_corrected <- (1-con)*p+(con*p_contaminant);
		d <- dmultinom(c(obs[,1], obs[,2]), prob=p_corrected, log=TRUE)
	}
	return(d)
}

## The function to calulate the likelihood over the entire set of parameters. It uses the previous function 'LikeCon'.
LikeFun <- function(observed, simulations, human_freq=NULL, contaminations=Contaminations, Binomial=FALSE) {
    require(parallel)
	nCores = detectCores()
	if (nCores < 1 | is.na(nCores)) { nCores=1 } # in case detectCores fails
	if(is.null(human_freq)) {
		L <- do.call("cbind", mclapply(contaminations, function(co) sapply(simulations, function(s) LikeCon(observed, sim=s[,c(1:2,5:6)], con=co, Binomial=Binomial)), mc.cores=nCores))
	} else {
		L <- do.call("cbind", mclapply(contaminations, function(co) sapply(simulations, function(s) LikeCon(observed, sim=cbind(s[,1:2], human_freq), con=co, Binomial=Binomial)), mc.cores=nCores))
	}
	## Likelihood values with contaminations in the columns and the evolutionary parameters (i.e. simulations) in the rows.
	return(L)
}

# Function to get the contamination value with the highest likelihood without a grid approach.
est_contam = function(obs, sim, hum, lower_bound=0, upper_bound=1) {
	# The default range for 'con' between 0 and 1
	sim = cbind(sim[rownames(obs),1:2], hum[rownames(obs),1:2])+1
	obs = obs[,1:2]+1
	likelihood <- function(con) {
		# Call the LikeCon function with the given 'con' value
		L <- LikeCon(obs, sim, con = con, Binomial = argv$Binomial)
		# Return the negative log-likelihood since optimize minimizes the function
		return(L)
	}
	# Perform optimization to find the value of 'con' that maximizes the likelihood
	result <- optimize(likelihood, interval = c(lower_bound, upper_bound), maximum = TRUE)
	# Extract the estimated value of 'con' with the highest likelihood
	return(c(c=result$maximum, logLike=result$objective))
}


# New function to calulate the likelihood and finding the maximum with optim. It uses the previous functions 'LikeCon' and 'est_contam'.
LikeFunOptim <- function(observed, simulations, human_freq=NULL, contam_low=0, contam_up=0) {
    require(parallel)
	nCores = detectCores()
	if (nCores < 1 | is.na(nCores)) { nCores=1 } # in case detectCores fails
	if(is.null(human_freq)) {
		# take human from simulations
		L <- do.call("rbind", mclapply(simulations, function(s) est_contam(observed[,1:2], sim=s[,1:2], hum=s[,5:6], lower_bound=contam_low, upper_bound=contam_up), mc.cores=nCores))
	} else {
		L <- do.call("rbind", mclapply(simulations, function(s) est_contam(observed[,1:2],  sim=s[,1:2], hum=human_freq[,1:2]), mc.cores=nCores))
	}
	# Likelihood and contamination in the columns and the evolutionary parameters (i.e. simulations) as row names.
	return(L)
}

# New function 
estimate <- function(mydata, expectations=Simulations, return_LLratio = FALSE, return_for_plot = FALSE, Binomial=FALSE) {
	# Function to estimate of split time and contamination.
	# The object 'Contaminations' is already present and 'Times' will be created from names(expectations). They will be the columns (Contaminations) and the rows (Times) of the 'res' matrix
	# If nothing is specified the source of contamination will be Yoruba.
	# ydata is the sum of all windows/bins together, which is used to calculate the likelihood. We thought of using the sum of the likelihoods of each window/bin but this turned out to not be the same when using the multinomial.
	ydata = mydata[[1]]
	# The parameters taken from the names of expectations (simulations)
	Times = t(sapply(strsplit(names(expectations), split="_"), as.numeric))
	colnames(Times) = c("T1", "T2", "pop")
	if(!exists("argv")) { # in case argv was not created.
		argv = NULL
	}
	if(!is.null(argv[["split_time_range"]])) {
		ids = which(Times[,"T1"] >= argv[["split_time_range"]][1] & Times[,"T1"] <= argv[["split_time_range"]][2])
		if (length(ids) > 2) {
			Times = Times[ids,]
			expectations = expectations[ids]
			rm(ids)
		} else {
			Warnings = c(Warnings, paste("(-s) split_time_range: No simulations within this range. Specify values between", paste(range(Times[,1]), collapse=" and "), "kya."))
		}
	}
	if(length(mydata) > 1) {
		ydata = Reduce('+', mydata)
	}
	# The old way:
	# res = LikeFun(observed=ydata[,1:2], simulations=expectations, human_freq=ydata[,3:4], contaminations=Contaminations, Binomial=Binomial)
	#max_LL = max(res); m1 <- which(res == max_LL, arr.ind=T); t1 <- Times[m1[,1],]; c1 <- Contaminations[m1[,2]]
	res = LikeFunOptim(observed=ydata[,1:2], simulations=expectations, human_freq=ydata[,3:4], contam_low=Contaminations[1], contam_up=Contaminations[2])
	m1 <- which.max(res[, "logLike"])
	max_LL = res[m1, "logLike"]
	t1 <- Times[m1,]
	c1 <- res[m1, "c"]
	if (return_LLratio) {
		m2 <- which(res[,"logLike"] >= max_LL-3.4)
		ll <- res[m2,"logLike"]
		t2 <- Times[m2,,drop=FALSE]
		c2 <- res[m2,"c"]
		out <- cbind(data.frame(round(t2)), Contamination=round(c2,4), logLikelihood=round(ll,4))[order(ll, decreasing=T),]
	} else {
		out <- c(round(t1), Contamination=round(c1,4), logLikelihood=round(max_LL,4))
		out <- data.frame(t(out))
	}
	if (return_for_plot) {
		return(list(Estimates=out, Results=res))
	} else {
		return(out)
	}
}

## Process the likehoods after running the function 'estimate'.
proc_est <- function(e_data, best_pop=NULL, quantile.prob = 0.95, comments = NULL) {
	# e_data must be a matrix sorted by the logLikelihood since the first will be taken as the point estimate.
	if (is.null(best_pop)) {
		best_pop=e_data[1,"pop"]
	}
	pop_lineages = list("vin"=c(1,5,6,7), "cha"=c(2,5,6,7), "alt"=c(3,6,7), "den"=c(4,7), "v-c"=c(1,2,5,6,7), nea=c(1:3,5:7), arc=c(1:7))
	nround = c(3,0,0,0,0,4,4,4,4)
	quantile.prob = (1-quantile.prob)/2
	quantile.prob = sort(c(quantile.prob, 1-quantile.prob))
	pop_best_lineage = pop_lineages[[best_pop]]
	if(nrow(e_data)>1) {
		tab = table(e_data[,"pop"])
		pops = names(tab)
		e_data = e_data[e_data[,"pop"] %in% pop_best_lineage, , drop=FALSE]
		est = apply(e_data[,-3],  2, function(x) quantile(x, prob=quantile.prob))
		out = data.frame(t(c(prop=nrow(e_data)/sum(tab), e_data[1,"T1"], est[,"T1"] , popSplit=e_data[1,"pop"], e_data[1,"Contamination"], est[,"Contamination"], e_data[1,"logLikelihood"])))
		colnames(out) = c("prop", "t","t_low", "t_high", "pop", "c", "c_low", "c_high", "logLike")
		if(sum(names(tab) %in% pop_best_lineage) < length(tab)) {
			comments = paste(paste(names(pop_lineages)[as.numeric(pops)],"=", round(tab/sum(tab)*100),"%",sep=""), collapse=",")
		}
	} else {
		out = data.frame(cbind(1,e_data[,c(1,1,1,3,4,4,4,5)]))
	}
	out <- data.frame(lapply(1:ncol(out), function(j) round(out[,j], nround[j])), stringsAsFactors = FALSE)
	colnames(out) = c("prop", "t","t_low", "t_high", "pop", "c", "c_low", "c_high", "logLike")
	return(list(Estimates=out, Comments=comments))
}


## Resample the table with rmultinom, assuming that all positions are indipendent.
resample <- function(tab, n=sum(tab)) {
	# n: total number of counts as size for the resampling.
	dm = dimnames(tab)
	if (ncol(tab) == 2) {
		tab[tab == 0] = 1e-7 # to avoid zero counts
		tab = matrix(rmultinom(1, n, c(tab[,1],tab[,2])/sum(tab)), ncol=2)
		dimnames(tab) = dimnames(dm)
	}
	return(tab)
}

#############################################################################
## Function(s) to generate the table given a certain type of ascertainment ##
#############################################################################

## Generate a table like the lineage assignment or class of sites from scrm simulations.
## This should replace the mktab functionfrom above.
sim2tab <- function(sim_data, human="africa", popN=c(europe=2,  africa=2,  vindija=2,  		chagyrskaya=2,  altai=2,  denisova=2,  test=2), contamination = 0, pop_cont_source="europe", pop_cont_sink="test", downsample=1, all_combinations = NULL, snp_array=NULL) {
	## the order of popN should not be changed.
	pop_a = names(popN)[3:6]
	pop_h = names(popN)[1:2]
	pop_ids = cumsum(popN)
	pop_ids = lapply(1:length(pop_ids), function(i) {x2=pop_ids[i]; x1=sum(as.numeric(pop_ids[i-1]))+1;x1:x2})
	names(pop_ids) = names(popN)
	if(is.null(all_combinations)) {
		all_combinations <- c("altai", "altai-chagyrskaya", "altai-chagyrskaya-denisova", "altai-denisova", "altai-vindija", "altai-vindija-denisova", "chagyrskaya", "chagyrskaya-denisova", "denisova", "human", "human-altai", "human-altai-chagyrskaya", "human-altai-chagyrskaya-denisova", "human-altai-denisova", "human-altai-vindija", "human-altai-vindija-denisova", "human-chagyrskaya", "human-chagyrskaya-denisova", "human-denisova", "human-neandertal", "human-vindija", "human-vindija-chagyrskaya", "human-vindija-chagyrskaya-denisova", "human-vindija-denisova", "neandertal", "neandertal-denisova", "vindija", "vindija-chagyrskaya", "vindija-chagyrskaya-denisova", "vindija-denisova")
	}
	# Traspose the data since using scrm within R was different than using scrm with the shell
	# This is done because it's faster to have many rows (snps) rather than many columns (haplotypes).
	sim_data <- lapply(sim_data, t)
    # Concatenate the simulations in order to match the real data.
    # Add contamination by replacing the sink population with the source in a given proportion of simulations. Notice that this is a cheap/fast way to do it because I mix simulations and not portion of sequences/individuals.
    if(contamination > 0) {
		n = length(sim_data)
        i <- sort(sample(1:n, round(n*as.numeric(contamination))))
        for (ii in i ) {
            x = sim_data[[ii]];
            x[, pop_ids[[pop_cont_sink]]] <- x[,pop_ids[[pop_cont_source]]]
            sim_data[[ii]] <- x
        }
    }
    sim_data = do.call("rbind", sim_data)
	# downsample the data if option downsample is a number between 0 and 1.
	if(is.numeric(downsample)) {
		if (downsample > 0 & downsample < 1) {
			n = nrow(sim_data)
			sim_data = sim_data[sample(1:n, round(n*downsample)), ]
		}
	} else {
		warning("downsample must be a number between 0 and 1.")
	}
	## Random sample one chromosome/haplotype for each population
    ids <- sapply(pop_ids[names(pop_ids) != human], function(i) if(length(i) > 1) {sample(i,1)} else {i})
    id_hum <- sample(pop_ids[[human]],2) # one to build the class of sites the other as a test
	## Substitute "human" with africa or europe depending which population was used for the lineage assignments.
    all_combinations <- sub("human", human, sub("neandertal", "altai-vindija-chagyrskaya", all_combinations))
    # Filter for the different ascertaimnets. It can still be improved and cleaned a bit.
    #snp_arrays=c("afr_high", "afr_unif", "afr_combined", "archaics_fixed", "archaics_var", "archaic_combined", "archaics_afr_combined")
	if(!is.null(snp_array)) {
		## With Africans
		daf_afr = rowSums(sim_data[,pop_ids[["africa"]]])/popN["africa"]
		f_bins=seq(0,1,0.02)
		f_ids = lapply(2:length(f_bins), function(i) which(daf_afr > f_bins[i-1] & daf_afr <= f_bins[i] ) )
		sfs_afr = sapply(f_ids,length)
		## Derive allele frequency in the archaics
		daf_arc = rowSums(sim_data[,unlist(pop_ids[pop_a]) ] )/sum(popN[pop_a])
		# Polymorphic in the Archaics
		arc_var = which(daf_arc > 0 & daf_arc < 1)
		# Fixed differences between archaics and Africans
		arc_fix1 = c(which(daf_arc == 0 & daf_afr == 1 ), which(daf_arc == 1 & daf_afr == 0 ))
		# Fixed derived in archaics, while not fixed derived in Africans
		arc_fix2 = which(daf_arc == 1 & daf_afr < 1)
		if (snp_array == "afr_high") { # SNPs having derived allele frequency in Africans > 10%
			snp_ids = which(daf_afr >= 0.1)
		}
		if (snp_array == "afr_unif") { # SNPs uniformly sampled in frequency-bins of 2%
			snp_ids = unlist(lapply(f_ids,function(x) sample(x,min(sfs_afr))))
		}
		if (snp_array == "afr_combined") { # Combination or better union of the above ones
			snp_ids = sort(unique(c(unlist(lapply(f_ids,function(x) sample(x,min(sfs_afr)))),  which(daf_afr >= 0.1))))
		}
		if (snp_array == "archaics_var") { # SNPs among the 4 archaic genomes
			snp_ids = arc_var
		}
		if (snp_array == "archaics_fixed") { # Fixed differences between Archaics and Africans (see above)
			snp_ids = sort(unique(c(afix1, afix2)))
		}
		if (snp_array == "archaics_combined") { # Union of the archaic snp_array
			snp_ids = sort(unique(c(arc_var, arc_fix1, arc_fix2)))
		}
		if (snp_array == "archaics_afr_combined") { # Union of all of them
			afr = c(unlist(lapply(f_ids, function(x) sample(x, min(sfs_afr)))),  which(daf_afr >= 0.1))
			snp_ids = sort(unique(c(afr, arc_var, arc_fix1, arc_fix2)))
		}
		if (snp_array == "archaics_plus" | snp_array == "archaic_plus") { # Newest array generated by Alba
			# Derived sites in one random African chromosome
			der1afr = which(sim_data[, sample(pop_ids[["africa"]],1)] ==1)
			# Fixed derived Archaics and not fixed derived in Africans
			arc_fix3 = which(daf_arc == 1 & daf_afr < 1)
			fix1  = which(daf_arc == 0 & daf_afr == 1)
			snp_ids = sort(unique(c(arc_var, fix1, arc_fix3, der1afr)))
		}
		sim_data = sim_data[snp_ids,]
		rm(snp_ids)
	}
	xdata = sim_data[,c(id_hum[1], ids[pop_a]) ]; colnames(xdata) = c(human, pop_a)
	xtest = sim_data[,c(id_hum[2], ids[!names(ids) %in% pop_a]) ]; colnames(xtest) = c("h1", "h2", "test")
	daf = rowSums(cbind(xdata, xtest))
	# Only select polymorphic sites
    polysites = which(daf > 0 & daf < ncol(xdata)+ncol(xtest))
	xdata = xdata[polysites,]
	xtest = xtest[polysites,]
	rm(sim_data, daf, polysites); for (i in 1:10) {g <- gc()}
    ## The informative sites ids for each of the classes given by "all_combinations" (an ugly code but efficient and fast)
    ids <- c()
    for (y in all_combinations) {
        i <- strsplit(y, split="-")[[1]]
        y1 = rowSums(xdata[, i, drop=F])
        y2 = rowSums(xdata[, !(colnames(xdata) %in% i), drop=F])
        ids[[y]] <- which(y1 == length(i) & y2 == 0)
    }
    ## get the count for the test population, and the two human populations.
    out = do.call("cbind", lapply(c("test", "h1", "h2"), function(p)  t(sapply(ids, function(i) if(length(i) > 0) { c(sum(xtest[i,p] == 0), sum(xtest[i,p] == 1))  } else { c(0,0) } )  )))
    return(out)
}

## As above but option random_sampling. So it can tabulate all types of genotypes if random_sampling is TRUE.
## There are also some other differences from above
sim2tab2 <- function(sim_data, human=names(popN)[1], popN=c(yoruba=2,  vindija=2,  chagyrskaya=2,  altai=2,  denisova=2,  test=1), outgroup=0, random_sampling=FALSE, contamination = 0, pop_cont_source=names(popN)[1], pop_cont_sink="test", pop_code=NULL, all_combinations = NULL, snp_array=NULL) {
	random_sampling = as.numeric(random_sampling)+1
	random_sampling[random_sampling> 2] = 2
	## the order of popN should not be changed.
	pops = names(popN)
	if(is.null(pop_code)) { # the last pop is consired always as a test.
		pop_code=rev(letters)[1:(length(pops)-1)]
	}
	pop_ids = cumsum(popN)
	pop_ids = lapply(1:length(pop_ids), function(i) {x2=pop_ids[i]; x1=sum(as.numeric(pop_ids[i-1]))+1;x1:x2})
	names(pop_ids) = names(popN)
	if(is.null(all_combinations)) {
		m = expand.grid(lapply(pop_code, function(x) c("", sapply(1:random_sampling, function(r) paste(x,r,sep="")))))
		all_combinations=apply(m,1,function(x) paste(x[x != ""],collapse="-"))
		all_combinations=all_combinations[all_combinations != ""]
	}
	## Traspose the data since using scrm within R was different than using scrm with the shell
	## This is done because it's faster to have many rows (snps) rather than many columns (haplotypes).
	sim_data <- lapply(sim_data, t)
    ## Concatenate the simulations in order to match the real data.
    ## Add contamination by replacing the sink population with the source in a given proportion of simulations. Notice that this is a cheap/fast way to do it because I mix simulations and not portion of sequences/individuals.
    if(contamination > 0) {
        i <- sort(sample(1:length(sim_data), round(length(sim_data)*as.numeric(contamination))))
        for (ii in i ) {
            x=sim_data[[ii]];
            x[, pop_ids[[pop_cont_sink]]] <- x[,pop_ids[[pop_cont_source]]]
            sim_data[[ii]] <- x
        }
    }
    sim_data = do.call("rbind", sim_data)
    ## Keep only sites where the outgroup is 0 (ancestral) and then remove it.
    if(length(outgroup) == 1 & outgroup > 0) {
		sim_data = sim_data[sim_data[,outgroup] == 0, -outgroup]
	}
    ## Filter for the different ascertaimnets. It can still be improved and cleaned a bit.
    #snp_arrays=c("afr_high", "afr_unif", "afr_combined", "archaics_fixed", "archaics_var", "archaic_combined", "archaics_afr_combined")
	if(!is.null(snp_array)) {
		## With Africans
		daf_afr = rowSums(sim_data[,pop_ids[["africa"]]])/popN["africa"]
		f_bins=seq(0,1,0.02)
		f_ids = lapply(2:length(f_bins), function(i) which(daf_afr > f_bins[i-1] & daf_afr <= f_bins[i] ) )
		sfs_afr = sapply(f_ids,length)
		## Derive allele frequency in the archaics
		daf_arc = rowSums(sim_data[,unlist(pop_ids[pop_a]) ] )/sum(popN[pop_a])
		# Polymorphic in the Archaics
		arc_var = which(daf_arc > 0 & daf_arc < 1)
		# Fixed differences between archaics and Africans
		arc_fix1 = c(which(daf_arc == 0 & daf_afr == 1 ), which(daf_arc == 1 & daf_afr == 0 ))
		# Fixed derived in archaics, while not fixed derived in Africans
		arc_fix2 = which(daf_arc == 1 & daf_afr < 1)
		if (snp_array == "afr_high") { # SNPs having derived allele frequency in Africans > 10%
			snp_ids = which(daf_afr >= 0.1)
		}
		if (snp_array == "afr_unif") { # SNPs uniformly sampled in frequency-bins of 2%
			snp_ids = unlist(lapply(f_ids,function(x) sample(x,min(sfs_afr))))
		}
		if (snp_array == "afr_combined") { # Combination or better union of the above ones
			snp_ids = sort(unique(c(unlist(lapply(f_ids,function(x) sample(x,min(sfs_afr)))),  which(daf_afr >= 0.1))))
		}
		if (snp_array == "archaics_var") { # SNPs among the 4 archaic genomes
			snp_ids = arc_var
		}
		if (snp_array == "archaics_fixed") { # Fixed differences between Archaics and Africans (see above)
			snp_ids = sort(unique(c(afix1, afix2)))
		}
		if (snp_array == "archaics_combined") { # Union of the archaic snp_array
			snp_ids = sort(unique(c(arc_var, arc_fix1, arc_fix2)))
		}
		if (snp_array == "archaics_afr_combined") { # Union of all of them
			afr = c(unlist(lapply(f_ids, function(x) sample(x, min(sfs_afr)))),  which(daf_afr >= 0.1))
			snp_ids = sort(unique(c(afr, arc_var, arc_fix1, arc_fix2)))
		}
		if (snp_array == "archaics_plus") { # Newest array generated by Alba
			# Fixed derived Archaics and not fixed derived in Africans
			arc_fix3 = which(daf_arc == 1 & daf_afr < 1)
			fix1  = which(daf_arc == 0 & daf_afr == 1)
			snp_ids = sort(unique(c(arc_var, fix1, arc_fix3)))
		}
		sim_data = sim_data[snp_ids,]
		rm(snp_ids)
	}
	## Random sample one chromosome/haplotype for each population
    ids <- lapply(pop_ids, function(i) if(length(i) > 1) {sample(i,random_sampling)} else {i})
	xtest = sim_data[,ids$test,drop=F];
    ids <- ids[1:(length(ids)-1)]
	xdata = do.call("cbind", lapply(ids, function(j) rowSums(sim_data[,j,drop=FALSE])))
	colnames(xdata) = pop_code
	colnames(xtest) = "test"
    ## The informative sites ids for each of the classes given by "all_combinations" (an ugly code but efficient and fast)
    out=c()
    for (y in all_combinations) {
        i <- strsplit(y, split="-")[[1]]
        ii = strsplit(i, split="")
        idx = sapply(ii,function(x) { z=as.numeric(x[2]); names(z)=x[1];z } )
		y1 =  rowSums(sapply(names(idx), function(x) xdata[,x] == idx[x]))
		y2 = rowSums(xdata[, !colnames(xdata) %in% (names(idx)), drop=F])
		ids = which(y1 == length(i) & y2 == 0)
        out[[y]] = unlist(lapply(colnames(xtest), function(i) c(Anc=sum(xtest[ids,i] == 0), Der=sum(xtest[ids,i] == 1) )))
    }
	if( random_sampling == 1) {
		all_combinations=gsub("1", "", all_combinations)
	}
    ## get the count for the test population, and the two human populations.
    out = do.call("rbind", out)
    rownames(out) = all_combinations
    return(out)
}


#########################################
## Function to estimate the parameters ##
## and to report the result in various formats to be processed later.

estimate_params <- function(mydata, pars=Parameters, expectations=Simulations, return_LLratio = FALSE, source_of_contamination = NULL, return_for_plot = FALSE) {
	## Function to estimate of split time and contamination.
	## The object 'Contaminations' is already present and 'Times' will be created from names(expectations). They will be the columns (Contaminations) and the rows (Times) of the res matrix
	# If nothing is specified the source of contamination will be Yoruba.
	# ydata is the sum of all windows/bins together, which is used to calculate the likelihood. We thought of using the sum of the likelihoods of each window/bin but this turned out to not be the same when using the multinomial.
	ydata = mydata[[1]]
	# First, remove parameters that do not vary. Actually it keeps the one that vary.
	pars=pars[, which(apply(pars,2,function(x) length(unique(x))) > 1)]
	# Filter/process the parameters
	if(!is.null(argv[["split_time_range"]])) {
		ids = which(pars[,"t_x"] >= argv[["split_time_range"]][1] & pars[,"t_x"] <= argv[["split_time_range"]][2])
		if (length(ids) > 2) {
			pars = pars[ids,]
			expectations = expectations[ids]
			rm(ids)
		} else {
			Warnings = c(Warnings, paste("(-s) split_time_range: No simulations within this range. Specify values between", paste(range(pars[,1]), collapse=" and "), "kya."))
		}
	}
	if(length(mydata) > 1) {
		for (y in 2:length(mydata)) { ydata = ydata + mydata[[y]] }
	}
	if(!is.null(argv$downsample)) {
		ydata = resample(ydata, argv$downsample)
	}
	PopInFreq = which(toupper(argv$population) == toupper(names(human_frequencies)))
	if(!is.null(argv$population) & length(PopInFreq) == 1) {
		source_of_contamination = human_frequencies[[PopInFreq]]
	} else {
		if(ncol(ydata) > 6) { # If there is the table of the SFS in the input use that as source  of contamination
			source_of_contamination = t(apply(ydata[,-(1:2)], 1, function(x) { d=sum(x*0:(length(x)-1)/(length(x)-1)); round(c(Ancestral=sum(x)-d, Derived=d)) }))
		}
	}
	res = LikeFun(observed=ydata[,1:2], simulations=expectations, human_freq=source_of_contamination, contaminations=Contaminations)
	## The ids where there is the highest likelihood: raw is for the split times (and population split) and column for contamination.
	max_LL = max(res)
	m1 <- which(res == max_LL, arr.ind=T)
	t1 <- pars[m1[,1],]
	c1 <- Contaminations[m1[,2]]
	if (nrow(m1) > 1) {
		# take the mean of the values when there is more than one maximum likelihood. Still I need to find a better solution, although this case is unlikely to happen.
		stop("More than one maximum likelihood.")
	}
	if (return_LLratio) {
		m2 <- which(res >= max_LL-3.4, arr.ind=T)
		ll <- res[m2]
		t2 <- pars[m2[,1],,drop=FALSE]
		c2 <- Contaminations[m2[,2]]
		out <- cbind(data.frame(t2), Contamination=round(c2,4),logLikelihood=round(ll,4))[order(ll, decreasing=T),]
	} else {
		out <- c(t1, Contamination=round(c1,4), logLikelihood=round(max_LL,4))
		out <- data.frame(t(out))
	}
	if (return_for_plot) {
		return(list(Estimates=out, Results=res))
	} else {
		return(out)
	}
}
