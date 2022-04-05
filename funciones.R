
## ---------------------------------------------------------------------------
## Teoria de portafolios - 2022-1
## Optimizacion de portafolios: Markowitz, Sharpe, Treynor
## y evaluación de desempeño
## Universidad Externado de Colombia
## Prof. Carlos Zapata [e-mail: cazapata25@gmail.com]
## @CpR- 2022 
## ---------------------------------------------------------------------------

## ---------------------------------------------------------------------------
## Importar precios Yahoo Finance
## Ej: para la informacion in-sample especificar:
## activos <- c("MA","AAPL","MSFT","GOOG","CVX","CAT","ABT","JNJ","MCD","HD") 
## fechai <- "2015-12-01"
## fechaf <- "2020-12-31"
## periodicidad <- "monthly" 

f.precios <- function(activos,fechai,fechaf,periodicidad){
    precios <- xts()
    for(i in 1:length(activos)){
        aux <- Ad(getSymbols(activos[i],from=fechai,to=fechaf,
                             periodicity=periodicidad,auto.assign=FALSE))
        aux <- na.approx(aux,na.rm=FALSE) # Interpolación de datos con NA
        precios <- cbind(precios,aux)
    }
    colnames(precios) <- activos
    tclass(precios) <- "Date"
    print("Los precios han sido importados correctamente")
    return(precios)
}

## ---------------------------------------------------------------------------
## ---------------------------------------------------------------------------
## Programa optimizacion Media-Varianza (MV)
## Teoria de portafolios - 2022-1
## Con cortos (con pesos negativos): short=1
## Sin cortos (sin pesos negativos): short=0
##-----------------------------------------------

modeloMV <- function(ret){
    # Inputs
    rf <- rf
    mu <- colMeans(ret)
    cov <- cov(ret)
    activos <- names(ret)
    # Optimizacion sin restricciones en cortos
    if(short == 1){
        ones <- rep(1,n)
        x <- t(mu)%*%solve(cov)%*%mu
        y <- t(mu)%*%solve(cov)%*%ones
        z <- t(ones)%*%solve(cov)%*%ones
        d <- x*z - y*y
        g <- (solve(cov,ones)%*%x-solve(cov,mu)%*%y)%*%solve(d)
        h <- (solve(cov,mu)%*%z-solve(cov,ones)%*%y)%*%solve(d)
        if(min(mu) > 0){rpmin = min(mu)*1.1}
        else{rpmin = 0.00}
        rpmax <- max(mu)*1.2
        nport <- 100
        j <- seq(rpmin,rpmax, length=nport) 
        wpo <- matrix(c(0), ncol=n, nrow=nport) 
        rpo <- matrix(c(0), nrow=nport)
        sigmapo <- matrix(c(0), nrow=nport)
        wj <- 0
        cont <- 1
        for(i in 1:nport){
            wj <- g + h*j[i] 
            wpo[cont,] <- t(wj)
            rpo[cont,] <- t(wj)%*%mu
            sigmapo[cont,] <- sqrt(t(wj)%*%cov%*%wj)
            cont <- cont+1
        }
        # PMVG
        wpmvg <- solve(cov,ones)%*%(1/z)
        rpmvg <- mu%*%wpmvg
        sigmapmvg <- sqrt(t(wpmvg)%*%cov%*%wpmvg)
        # Sharpe
        Er <- mu-rf 
        Z <- solve(cov,Er)  
        sumZ <- sum(Z) 
        wpt <- Z/sumZ 
        rpt <- t(wpt)%*%mu
        sigmapt <- sqrt(t(wpt)%*%cov%*%wpt)
        wpmvg <- t(wpmvg)
        wpt <- t(wpt)
        
        MV <- list()
        MV[[1]] <- wpo
        MV[[2]] <- rpo
        MV[[3]] <- sigmapo
        MV[[4]] <- t(wpmvg)
        MV[[5]] <- rpmvg
        MV[[6]] <- sigmapmvg
        MV[[7]] <- t(wpt)
        MV[[8]] <- rpt 
        MV[[9]] <- sigmapt
        return(MV)
    }
    # Con restricciones en corto
    else {
        # FE    
        library(quadprog)
        if(min(mu) > 0){rpmin = min(mu)*1.1}
        else{rpmin = 0.00}
        rpmax <- max(mu)*0.99
        n <- length(mu)
        nport <- 1000
        j <- seq(rpmin,rpmax,length=nport)
        sigmapo <- matrix(0,nrow=nport)
        wpo <- matrix(0,nrow=nport, ncol=n)
        Amat <- t(rbind(rep(1,n),mu,diag(1,nrow=n)))
        dvec <- rep(0,n) 
        Dmat <- 2*cov
        for(i in 1:nport){
            bvec <- c(1,j[i],rep(0,n))
            result <- solve.QP(Dmat=Dmat,dvec=dvec,Amat=Amat,bvec=bvec,meq=2)
            wpo[i,] <- result$solution
            sigmapo[i,] <- sqrt(result$value)
        }
        rpo <- j
        colnames(wpo) <- c(activos)
        # PMVG
        pmvg <- cbind(sigmapo,wpo)
        pmvg.sort <- pmvg[order(pmvg[,1]),]
        pmvg.sel <- cbind(pmvg.sort[1,])
        wpmvg <- cbind(round(pmvg.sel[2:length(pmvg.sel)],6))
        rownames(wpmvg) <- c(activos)
        rpmvg <- mu%*%wpmvg
        sigmapmvg <- sqrt(t(wpmvg)%*%cov%*%wpmvg)
        # Sharpe    
        sharpe_port <- (rpo-rf)/sigmapo
        sharpe <- cbind(sharpe_port,wpo)
        sharpe.sort <- sharpe[order(-sharpe[,1]),]
        sharpe.sel <- cbind(sharpe.sort[1,])
        wpt <- round(cbind(sharpe.sel[2:length(sharpe.sel)]),6)
        rownames(wpt) <- c(activos)
        rpt <- mu%*%wpt
        sigmapt <- sqrt(t(wpt)%*%cov%*%wpt)
        
        MV <- list()
        MV[[1]] <- wpo
        MV[[2]] <- rpo
        MV[[3]] <- sigmapo
        MV[[4]] <- wpmvg
        MV[[5]] <- rpmvg
        MV[[6]] <- sigmapmvg
        MV[[7]] <- wpt
        MV[[8]] <- rpt 
        MV[[9]] <- sigmapt
        return(MV)
    }
}

## ---------------------------------------------------------------------------
## ---------------------------------------------------------------------------
## Programa optimización de Treynor
## -----------------------------------
## Recomendacion:usar minimo 5 años
## de historia en periodicidad mensual
## -----------------------------------

m.treynor <- function(retornos,r.indice){
    n <- ncol(retornos)
    betas <- matrix(0,ncol=n)
    varerror<- matrix(0,ncol=n)
    # Regresion iterativa para los parametros
    for(i in 1:n){
        modelo <- lm(retornos[,i]~r.indice)
        betas[i] <- modelo[["coefficients"]][2]
        varerror[i] <- var(modelo[["residuals"]])
    }
    treynori <- (mu-rf)/betas
    # Calculo de los ratios 1 y 2 y las sumas acumuladas
    matriz <- t(rbind(treynori,betas,varerror,mu,sigma))
    colnames(matriz) <- c("Treynor","Betas","VaError","Mu","Sigma")
    matriz.ord <- matriz[order(-matriz[,1]),]
    ratio1 <- ((matriz.ord[,4]-rf)*matriz.ord[,2])/matriz.ord[,3]
    ratio2 <- matriz.ord[,2]^2/matriz.ord[,3]
    suma1 <- cumsum(ratio1)
    suma2 <- cumsum(ratio2)
    sigmam <- sd(r.indice)
    tasac <- (sigmam^2*suma1)/(1+sigmam^2*suma2)
    diff <- matriz.ord[,1] - tasac
    cond.diff <- diff[ diff>0 ]
    n.optimo <- length(cond.diff)
    cmax <- max(tasac)
    zi <- (matriz.ord[,2]/matriz.ord[,3])*(matriz.ord[,1]-cmax)
    zi <- pmax(zi,0)
    wpot <- zi/sum(zi)
    wpot <- rbind(wpot[activos])
    wpot <- t(wpot)
    rpot <- t(wpot)%*%mu
    sigmapot <- sqrt(t(wpot)%*%cov%*%wpot)
    
    MT <- list()
    MT[[1]] <- wpot
    MT[[2]] <- rpot
    MT[[3]] <- sigmapot
    MT[[4]] <- n.optimo
    return(MT)
}


##-----------------------------------------------
## Programa optimizacion Semi Media-Varianza (SMV)
## Teoria de portafolios - 2022-1
## Con cortos (con pesos negativos): short=1
## Sin cortos (sin pesos negativos): short=0
##-----------------------------------------------

m.sortino <- function(retornos,h){
    # Inputs
    h <- h
    rf <- rf
    mu <- colMeans(retornos)
    semiretornos <- pmin(retornos,h)
    cov <- cov(semiretornos) # Semi-covarianzas
    activos <- names(retornos)

    # Optimizacion sin restricciones en cortos
    if(short == 1){
        ones <- rep(1,n)
        x <- t(mu)%*%solve(cov)%*%mu
        y <- t(mu)%*%solve(cov)%*%ones
        z <- t(ones)%*%solve(cov)%*%ones
        d <- x*z - y*y
        g <- (solve(cov,ones)%*%x-solve(cov,mu)%*%y)%*%solve(d)
        h <- (solve(cov,mu)%*%z-solve(cov,ones)%*%y)%*%solve(d)
        if(min(mu) > 0){rpmin = min(mu)*1.1}
        else{rpmin = 0.00}
        rpmax <- max(mu)*1.2
        nport <- 100
        j <- seq(rpmin,rpmax, length=nport) 
        wpo <- matrix(c(0), ncol=n, nrow=nport) 
        rpo <- matrix(c(0), nrow=nport)
        sigmapo <- matrix(c(0), nrow=nport)
        wj <- 0
        cont <- 1
        for(i in 1:nport){
            wj <- g + h*j[i] 
            wpo[cont,] <- t(wj)
            rpo[cont,] <- t(wj)%*%mu
            sigmapo[cont,] <- sqrt(t(wj)%*%cov%*%wj)
            cont <- cont+1
        }
        # PMVG
        wpmvg <- solve(cov,ones)%*%(1/z)
        rpmvg <- mu%*%wpmvg
        sigmapmvg <- sqrt(t(wpmvg)%*%cov%*%wpmvg)
        # Sortino
        Er <- mu-rf 
        Z <- solve(cov,Er)  
        sumZ <- sum(Z) 
        wpt <- Z/sumZ 
        rpt <- t(wpt)%*%mu
        sigmapt <- sqrt(t(wpt)%*%cov%*%wpt)
        wpmvg <- t(wpmvg)
        wpt <- t(wpt)
        
        SMV <- list()
        SMV[[1]] <- wpo
        SMV[[2]] <- rpo
        SMV[[3]] <- sigmapo
        SMV[[4]] <- t(wpmvg)
        SMV[[5]] <- rpmvg
        SMV[[6]] <- sigmapmvg
        # Punto tangente de Sortino 7-9
        SMV[[7]] <- t(wpt)
        SMV[[8]] <- rpt 
        SMV[[9]] <- sigmapt
        return(SMV)
    }
    # Con restricciones en corto
    else {
    # FE    
    library(quadprog)
    if(min(mu) > 0){rpmin = min(mu)*1.1}
    else{rpmin = 0.00}
    rpmax <- max(mu)*0.99
    n <- length(mu)
    nport <- 1000
    j <- seq(rpmin,rpmax,length=nport)
    sigmapo <- matrix(0,nrow=nport)
    wpo <- matrix(0,nrow=nport, ncol=n)
    Amat <- t(rbind(rep(1,n),mu,diag(1,nrow=n)))
    dvec <- rep(0,n) 
    Dmat <- 2*cov
    for(i in 1:nport){
      bvec <- c(1,j[i],rep(0,n))
      result <- solve.QP(Dmat=Dmat,dvec=dvec,Amat=Amat,bvec=bvec,meq=2)
      wpo[i,] <- result$solution
      sigmapo[i,] <- sqrt(result$value)
    }
    rpo <- j
    colnames(wpo) <- c(activos)
    # PMVG
    pmvg <- cbind(sigmapo,wpo)
    pmvg.sort <- pmvg[order(pmvg[,1]),]
    pmvg.sel <- cbind(pmvg.sort[1,])
    wpmvg <- cbind(round(pmvg.sel[2:length(pmvg.sel)],6))
    rownames(wpmvg) <- c(activos)
    rpmvg <- mu%*%wpmvg
    sigmapmvg <- sqrt(t(wpmvg)%*%cov%*%wpmvg)
    # Sortino
    sortino_port <- (rpo-rf)/sigmapo
    sortino <- cbind(sortino_port,wpo)
    sortino.sort <- sortino[order(-sortino[,1]),]
    sortino.sel <- cbind(sortino.sort[1,])
    wpt <- round(cbind(sortino.sel[2:length(sortino.sel)]),6)
    rownames(wpt) <- c(activos)
    rpt <- mu%*%wpt
    sigmapt <- sqrt(t(wpt)%*%cov%*%wpt)

    SMV <- list()
    SMV[[1]] <- wpo
    SMV[[2]] <- rpo
    SMV[[3]] <- sigmapo
    SMV[[4]] <- wpmvg
    SMV[[5]] <- rpmvg
    SMV[[6]] <- sigmapmvg
    SMV[[7]] <- wpt
    SMV[[8]] <- rpt 
    SMV[[9]] <- sigmapt
    return(SMV)
    }
}


## -----------------------------------
## Funcion de optimizacion medida Omega
## -----------------------------------

m.omega <- function(retornos,h){
    library("ROI")
    library("ROML")
    library("ROML.portfolio")
    h <- h
    short <- short
    if(short == 1){lb=-1}
    else{lb=0}
    m <- model()
    m$variable(portfolio, lb = lb) # the portfolio choice vector; 
    m$maximize( omega(portfolio) )
    opt <- optimize(m, solver="glpk", 
                    data=list(returns = coredata(retornos))) 
    wpomega <- round(opt$solution[grep("portfolio", names(opt$solution))]/
                         opt$solution[grep("z", names(opt$solution))], 6)
    
    names(wpomega) <- activos
    rpomega <- mu%*%wpomega
    sigmapomega <- sqrt(t(wpomega)%*%cov%*%wpomega)
    
    PO <- list()
    PO[[1]] <- cbind(wpomega)
    PO[[2]] <- rpomega
    PO[[3]] <- sigmapomega
    return(PO)
}

## -----------------------------------
## FUncion de optimizacion del CVaR
## -----------------------------------

m.cvar <- function(retornos,alpha){
    library("ROI")
    library("ROML")
    library("ROML.portfolio")
    alpha <- alpha
    short <- short
    if(short==1){lb=-1}
    else{lb=0}
    
    m <- model()
    m$variable(portfolio, lb = lb) # the portfolio choice vector; 
    m$minimize( cvar(portfolio, alpha) )
    m$subject_to( budget_norm(portfolio) )
    opt <- optimize(m, solver="glpk", 
                    data=list(returns = coredata(retornos)))

    wpcvar <- round(opt$solution[grep("portfolio", names(opt$solution))], 4)
    names(wpcvar) <- activos
    rpcvar <- mu%*%wpcvar
    sigmapcvar <- sqrt(t(wpcvar)%*%cov%*%wpcvar)
    
    MCVAR <- list()
    MCVAR[[1]] <- wpcvar
    MCVAR[[2]] <- rpcvar
    MCVAR[[3]] <- sigmapcvar
    return(MCVAR)
}

# ---------------------------------------
# Risk Parity optimization
# Forma 1: solucion del polinomio
# ------------------------------------------

m.riskparity <- function(cov){
    library(rootSolve)
    b <- rep(1/n,n)
    # Se cosntruye el polinomio para el ejercicio iterativo
    f_root <- function(x,coef){
        cov <- coef
        n <- ncol(cov)
        return(cov%*%x - b/x)
    }
    # Se obtienen las raices del polinomio
    x_root <- multiroot(f_root, start=b, parms=cov)$root
    wRP <- x_root/sum(x_root)
    names(wRP) <- activos
    
    rpRP <- mu%*%wRP
    sigmapRP <- sqrt(t(wRP)%*%cov%*%wRP)
    CR <- wRP*(cov%*%wRP)
    CR <- round(CR,6)
    
    RPP <- list()
    RPP[[1]] <- cbind(wRP)
    RPP[[2]] <- rpRP
    RPP[[3]] <- sigmapRP
    RPP[[4]] <- CR
    return(RPP)
}


## ---------------------------------------------------------------------------
## ---------------------------------------------------------------------------
## Programa Performance Evaluation
## -----------------------------------------------
## Recomendacion input: vector de pesos en columna

performance <- function(retornos,pesos){
    t <- nrow(retornos)
    n <- ncol(retornos)
    rport <- matrix(0,nrow=t,ncol=1)
    colnames(rport) <- c("Retornos")
    vport <- matrix(0,nrow=t,ncol=1)
    colnames(vport) <- c("Valor")
    pesos <- pesos 
    
    # Retornos historicos
    ret.hist <- retornos%*%pesos
    rport[,1] <- ret.hist
    
    # Valor del portafolio
    port.v <- matrix(0, nrow=t)
    port.v[1] <- valor
    for(i in 2:t){
        port.v[i] <- port.v[i-1]*exp(ret.hist[i-1])
    }
    vport[,1] <- port.v
    
    DH <- list()
    DH[[1]] <- vport
    DH[[2]] <- rport
    return(DH)
}


## Cont...


## End.
## ---------------------------------------------------------------------------



