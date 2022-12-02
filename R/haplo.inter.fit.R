`haplo.inter.fit` <-
function(geno, var2, dep, adj = NULL, fam, haplo.freq.min, ...)
{
	param <- "var2"; #This is the name of the formal parameter that contains the interaction variable

	# haplo.inter(genotip, Datos$sexo, Datos$grupo)
	# haplo.inter(genotip, Datos$sexo, Datos$grupo, Datos[,c("rcal.dia", "n.edad")],0.005,1)

	# Funcio que retorna les tres taules d'interaccions a partir dels coef i se de l'ajust d'un model haplo.glm

	# Li passo els seguents par?metres per crear un model haplo.glm amb interaccions

	# La variable geno es del tipus: (genotip <- setupGeno(cbind(allele(g.81,1:2), allele(g.82,1:2), allele(g.83,1:2))) on g.81<-genotype(XRCC1.81)  )
	# var2 (Variable que interacciona amb l'haplotip)
	# dep (Variable depenent)
	# adj (Variables d'ajust)
	# fam (funci? pel glm "binomial"/"gaussian")
	# haplo.freq.min (haplo.freq.min = 0.01) Frequ?ncia m?nima haplot?pica.
	# num.status (variable resposta categ?rica o num?rica?)

	# Models
	if (is.null(adj))
	{
		mod <- haplo.stats::haplo.glm( dep~geno*var2, family=fam, allele.lev=attributes(geno)$unique.alleles, 
		                               control = haplo.stats::haplo.glm.control(haplo.freq.min=haplo.freq.min,...))
		mod.b <- haplo.stats::haplo.glm( dep~geno+var2, family=fam, allele.lev=attributes(geno)$unique.alleles, 
		                                 control = haplo.stats::haplo.glm.control(haplo.freq.min=haplo.freq.min,...))
	}
	else
	{
		mod <- haplo.stats::haplo.glm( dep~.+geno*var2, family=fam, allele.lev = attributes(geno)$unique.alleles, 
		                               control = haplo.stats::haplo.glm.control(haplo.freq.min=haplo.freq.min,...), data=adj)
		mod.b <- haplo.stats::haplo.glm( dep~.+geno+var2, family=fam, allele.lev=attributes(geno)$unique.alleles, 
		                                 control = haplo.stats::haplo.glm.control(haplo.freq.min=haplo.freq.min,...), data=adj)
	}

	# Calculo la p d'interacci? entre l'haplotip i la variable
	pval.haplo <- 1-pchisq(mod$lrt$lrt-mod.b$lrt$lrt, mod$lrt$df-mod.b$lrt$df)

	#  Guardo els coeficients en el vector co
	co <- mod$coef
	noms.coef <- names(co)

	# Matriu de covari?ncies
	# Selecciono les mateixes files i columnes que a la matriu de coef, ja que la matriu q retorna mod$var.mat no t? dimnames i estan tots els haplotips sense agrupar els rare. 

	mat.cov <- mod$var.mat[1:length(mod$coef), 1:length(mod$coef)]
	dimnames(mat.cov)<- list(names(mod$coef), names(mod$coef))

	# Creo la nova matriu (mat.coef) amb els resultats dels coeficients, per files var2 i columnes geno, amb el format de taula final.
	mat.coef <- matrix(nrow=1+length(mod$haplo.names), ncol=length(levels(var2)))
	dimnames(mat.coef) <- list(c(mod$haplo.base, mod$haplo.names), levels(var2))

	# Creo una matriu amb les vari?ncies del mateix format que la de coeficients.
	mat.var <- mat.coef
	mat.coef[1,1] <- 0 # Refer?ncia

	# creo matrius per taules marginals de row i col
	mat.coef.c <- mat.coef
	mat.var.c <- mat.var
	mat.coef.r <- mat.coef
	mat.var.r <- mat.var
 
	# Inicialitzo i, j
	i <- 1
	j <- 1

	# Selecciono la posici? dels efectes principals (de tot el conjunt elimino les interaccions)
	# Omplo la primera columna (correspon als efectes principals de la geno)
	for (i in 1:(dim(mat.coef)[1]-1))
	{
		pos <- rownames(mat.coef)[i+1];

		mat.coef[i+1,1] <- co[pos];
		mat.var[i+1,1] <- mat.cov[pos,pos];
		mat.coef.c[i+1,1] <- co[pos]
		mat.var.c[i+1,1] <- mat.cov[pos,pos]
		mat.coef.r[i+1,1] <- NA
		mat.var.r[i+1,1] <- NA
	}

	# Omplo la primera fila (correspon als efectes principals de la var2)
	for (j in 1:(dim(mat.coef)[2]-1))
	{
		pos <- paste(param,colnames(mat.coef)[j+1],sep="");

		mat.coef[1,j+1] <- co[pos];
		mat.var[1,j+1] <- mat.cov[pos,pos];

		mat.coef.c[1,j+1] <- NA
		mat.var.c[1,j+1] <- NA

		mat.coef.r[1,j+1] <- co[pos]
		mat.var.r[1,j+1] <- mat.cov[pos,pos]
	}

	# Omplo les interaccions

	for (j in 1:(dim(mat.coef)[2]-1))
	{
		for (i in 1:(dim(mat.coef)[1]-1))
		{

			pos.x <- rownames(mat.coef)[i+1];
			pos.y <- paste(param,colnames(mat.coef)[j+1],sep="");
			pos.xy <- paste(rownames(mat.coef)[i+1],":",param,colnames(mat.coef)[j+1],sep="");

			mat.coef[i+1,j+1] <- mat.coef[i+1,1] + mat.coef[1,j+1] + co[pos.xy]
			mat.var[i+1,j+1]  <- mat.var[i+1,1] + mat.var[1,j+1] + mat.cov[pos.xy, pos.xy] + 2*(mat.cov[pos.x, pos.y] + mat.cov[pos.x, pos.xy] + mat.cov[pos.y, pos.xy])

			mat.coef.c[i+1,j+1] <- mat.coef.c[i+1,1] + co[pos.xy]
			mat.var.c[i+1,j+1]  <- mat.var.c[i+1,1] +  mat.cov[pos.xy, pos.xy] + 2*mat.cov[pos.x, pos.xy]

			mat.coef.r[i+1,j+1] <- mat.coef.r[1,j+1] + co[pos.xy]
			mat.var.r[i+1,j+1]  <- mat.var.r[1,j+1] + mat.cov[pos.xy, pos.xy] + 2*mat.cov[pos.y, pos.xy]
		}
	}

	# Calculo la matriu d'ORs/coefs i ICs
	mat.or <- mat.coef
	mat.li <- mat.coef - (1.96 * sqrt(mat.var))
	mat.ls <- mat.coef + (1.96 * sqrt(mat.var))

	mat.or.c <- mat.coef.c
	mat.li.c <- mat.coef.c - (1.96 * sqrt(mat.var.c))
	mat.ls.c <- mat.coef.c + (1.96 * sqrt(mat.var.c))

	mat.or.r <- mat.coef.r
	mat.li.r <- mat.coef.r - (1.96 * sqrt(mat.var.r))
	mat.ls.r <- mat.coef.r + (1.96 * sqrt(mat.var.r))

	#En cas de regressi? log?stica transformem els coeficients a OR's
	if (fam=="binomial")
	{
		mat.or <- exp(mat.or)
		mat.li <- exp(mat.li)
		mat.ls <- exp(mat.ls)

		mat.or.c <- exp(mat.or.c)
		mat.li.c <- exp(mat.li.c)
		mat.ls.c <- exp(mat.ls.c)

		mat.or.r <- exp(mat.coef.r)
		mat.li.r <- exp(mat.coef.r - 1.96 * sqrt(mat.var.r))
		mat.ls.r <- exp(mat.coef.r + 1.96 * sqrt(mat.var.r))
	}

	# Ordeno les tres matrius d'ORs i IC en una taula final (mat.fi)
	mat.fi <- NULL
	mat.fi.c <- NULL
	mat.fi.r <- NULL
	i <- 1
	j <- 1

	for(i in 1:length(levels(var2)))
	{
		mat.fi <- cbind(mat.fi, mat.or[,i], mat.li[,i], mat.ls[,i])
		dimnames(mat.fi)[[2]][j:(j+2)] <- c(dimnames(mat.or)[[2]][i], "li", "ls")

		mat.fi.c <- cbind(mat.fi.c, mat.or.c[,i], mat.li.c[,i], mat.ls.c[,i])
		dimnames(mat.fi.c)[[2]][j:(j+2)] <- c(dimnames(mat.or.c)[[2]][i], "li", "ls")
		mat.fi.r <- cbind(mat.fi.r, mat.or.r[,i], mat.li.r[,i], mat.ls.r[,i])
		dimnames(mat.fi.r)[[2]][j:(j+2)] <- c(dimnames(mat.or.r)[[2]][i], "li", "ls")

		j <- j + 3
	}

	# Afegeixo una primera columna a mat.fi amb la frequ?ncia dels haplotips
	HapFreq <- NA
	mat.fi <- cbind(HapFreq, mat.fi)
	mat.fi.c <- cbind(HapFreq, mat.fi.c)
	mat.fi.r <- cbind(HapFreq, mat.fi.r)

	# Arreglo les etiquetes de les files de la matriu corresponent als haplotips.
	# Elimino de les etiquetes el "geno.", de forma que em quedi nom?s el n?mero (geno.1 per 1),
	# per aix? selecciono a partir de la posici? 6 i fins la 9 pensant que el m?xim de llargada s?n els rare.
	# I selecciono a partir de la segona fila pq a la primera hi ha el base que ja ?s sempre un n?mero.

	dimnames(mat.fi)[[1]][2:nrow(mat.fi)] <- substr(dimnames(mat.fi)[[1]][2:nrow(mat.fi)], 6,9)
	dimnames(mat.fi.c)[[1]][2:nrow(mat.fi.c)] <- substr(dimnames(mat.fi.c)[[1]][2:nrow(mat.fi.c)], 6,9)
	dimnames(mat.fi.r)[[1]][2:nrow(mat.fi.r)] <- substr(dimnames(mat.fi.r)[[1]][2:nrow(mat.fi.r)], 6,9)

	# Un cop tinc els dimnames amb n?meros (3,4,...), li ajunto les etiquetes i freq. haplot?piques que li corresponen.
	# Ex: geno.2, ara ?s 2, i correspon a CGA, freq 0.34.
	i <- 1

	for(i in 1:nrow(mat.fi))
	{
		if (dimnames(mat.fi)[[1]][i]!="rare")
		{
			num <- as.numeric(dimnames(mat.fi)[[1]][i]) 
			mat.fi[i,c("HapFreq")] <- mod$haplo.freq[num]
			dimnames(mat.fi)[[1]][i] <- paste(mod$haplo.unique[num,], collapse="")

			mat.fi.c[i,c("HapFreq")] <- mod$haplo.freq[num]
			dimnames(mat.fi.c)[[1]][i] <- paste(mod$haplo.unique[num,], collapse="")
			mat.fi.r[i,c("HapFreq")] <- mod$haplo.freq[num]
			dimnames(mat.fi.r)[[1]][i] <- paste(mod$haplo.unique[num,], collapse="")
		}
 		else if (dimnames(mat.fi)[[1]][i]=="rare")
 		{
			mat.fi[i,c("HapFreq")] <- 1-sum(mat.fi[1:(nrow(mat.fi)-1),1])
			mat.fi.c[i,c("HapFreq")] <- 1-sum(mat.fi.c[1:(nrow(mat.fi.c)-1),1])
			mat.fi.r[i,c("HapFreq")] <- 1-sum(mat.fi.r[1:(nrow(mat.fi.r)-1),1])
		}
	}

	# Ordeno la taula final per frequ?ncia haplot?pica (de major a menor), deixant els rare sempre al final.
	# I suposant que el haplobase sempre ser? el m?s frequent i quedar? a la primera categoria ( en cas que no fos aix? no passaria res, nom?s que 
	# la categoria de refer?ncia no es mostraria a la primera fila.

	mat.fi <- mat.fi[c(order(mat.fi[1:(nrow(mat.fi)-1),c("HapFreq")], decreasing=TRUE), nrow(mat.fi)),]
	mat.fi.c <- mat.fi.c[c(order(mat.fi.c[1:(nrow(mat.fi.c)-1),c("HapFreq")], decreasing=TRUE), nrow(mat.fi.c)),]
	mat.fi.r <- mat.fi.r[c(order(mat.fi.r[1:(nrow(mat.fi.r)-1),c("HapFreq")], decreasing=TRUE), nrow(mat.fi.r)),]

	if (fam=="binomial")
	{
		# Si hi ha algun valor molt gran a la taula final li assigno el valor Inf.
		mat.fi[mat.fi>999] <- Inf
		mat.fi.c[mat.fi.c>999] <- Inf
		mat.fi.r[mat.fi.r>999] <- Inf
	}

	# Arrodoneixo a 4 decimals la frequ?ncia haplot?pica i la resta a 2.
	mat.fi[, c("HapFreq")] <- round(mat.fi[, c("HapFreq")], 4)
	mat.fi[, 2:ncol(mat.fi)] <- round(mat.fi[, 2:ncol(mat.fi)], 2)

	mat.fi.c[, c("HapFreq")] <- round(mat.fi.c[, c("HapFreq")], 4)
	mat.fi.c[, 2:ncol(mat.fi.c)] <- round(mat.fi.c[, 2:ncol(mat.fi.c)], 2)

	mat.fi.r[, c("HapFreq")] <- round(mat.fi.r[, c("HapFreq")], 4)
	mat.fi.r[, 2:ncol(mat.fi.r)] <- round(mat.fi.r[, 2:ncol(mat.fi.r)], 2)

	list(mat.fi=mat.fi, mat.fi.c=mat.fi.c, mat.fi.r=mat.fi.r, pval=pval.haplo) 
	
}

