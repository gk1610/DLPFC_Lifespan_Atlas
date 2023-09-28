library(parallel)
library(dplyr)

#the least square cosinor package--very long ends at line 433
one_cosinor_OLS = function(tod = truth$tod, y = exp.mat[502, ], alpha = 0.05){
  #alpha is the critical level for equal tailed CI
  n = length(tod)
  x1 = cos(2*pi*tod/24)
  x2 = sin(2*pi*tod/24)
  
  
  # mat.X = matrix(c(rep(1, n), x1, x2), ncol = 3, byrow = FALSE)
  # mat.XX = t(mat.X)%*%mat.X#mat.XX = mat.S
  mat.S = matrix(c(n, sum(x1), sum(x2), 
                   sum(x1), sum(x1^2), sum(x1*x2), 
                   sum(x2), sum(x1*x2), sum(x2^2)), 
                 nrow = 3, byrow = TRUE)
vec.d = c(sum(y), sum(y*x1), sum(y*x2))

  
  mat.S.inv = solve(mat.S)  
  est = mat.S.inv%*%vec.d
  m.hat = est[1]
  beta1.hat = est[2]
  beta2.hat = est[3]
 # truth$phase[2]; truth$amplitude[2]; truth$M[502]
  A.hat = sqrt(beta1.hat^2 + beta2.hat^2)
  
  phase.res = get_phase(beta1.hat, beta2.hat)
  phase.hat = phase.res$phase

  #inference 
  TSS = sum((y-mean(y))^2)
  yhat = m.hat + beta1.hat*x1+beta2.hat*x2
  RSS = sum((y-yhat)^2)
  MSS = TSS-RSS
  Fstat = (MSS/2)/(RSS/(n-3))
  pval = pf(Fstat, 2, n-3, lower.tail = FALSE)
  
  #CI (M and se for A and phi)
  sigma2.hat = RSS/(n-3)
  sigma.hat = sqrt(sigma2.hat)
  CI.m.hat.radius = qt(1-alpha/2, n-3)*sigma.hat*mat.S.inv[1, 1]
  se.A.hat = sigma.hat*sqrt(mat.S.inv[2, 2]^cos(phase.hat)^2
                            -2*mat.S.inv[2, 3]*sin(phase.hat)*cos(phase.hat)
                            +mat.S.inv[3, 3]*sin(phase.hat)^2)
  se.phase.hat = sigma.hat*sqrt(mat.S.inv[2, 2]^sin(phase.hat)^2
                              +2*mat.S.inv[2, 3]*sin(phase.hat)*cos(phase.hat)
                              +mat.S.inv[3, 3]*cos(phase.hat)^2)/A.hat
  
  #CI (derive conservative CI for phi)
  B11 = sum((x1-mean(x1))^2)
  B12 = sum((x1-mean(x1))*(x2-mean(x2)))
  B22 = sum((x2-mean(x2))^2)
  C1 = -(B11*beta1.hat+B12*beta2.hat)/(B22*beta2.hat+B12*beta1.hat)
  C2 = -(2*sigma2.hat*qf(1-alpha, 2, n-3)-2*B12*beta1.hat*beta2.hat-B11*beta1.hat^2-B22*beta2.hat^2)/(B22*beta2.hat+B12*beta1.hat)
  D1 = B22*C1^2+B11+2*B12*C1
  D2 = 2*B22*C1*C2+2*B12*C2-(B12*beta1.hat+B22*beta2.hat)*C1-B11*beta1.hat-B12*beta2.hat
  D3 = B22*C2^2-(B12*beta1.hat+B22*beta2.hat)*C2
  
  #calculate CI of phi
  #check if 0 is in ellipse 
  zero.in.ellipse = B11*beta1.hat^2+2*B12*beta1.hat*beta2.hat+B22*beta2.hat^2 < 2*sigma2.hat*qf(1-alpha, 2, n-3)
  if(!zero.in.ellipse){
    delta.poly = D2^2-4*D1*D3
    if(delta.poly<0){
      phi.lower.limit = list(tan = -99, phase = -99)
      phi.upper.limit = list(tan = -99, phase = -99)
    }else{
      phi.beta1.roots = c((-D2-sqrt(delta.poly))/(2*D1), (-D2+sqrt(delta.poly))/(2*D1))
      phi.beta2.roots = C1*phi.beta1.roots+C2
    
      # phi.limit1 = get_phase(phi.beta1.roots[1], phi.beta2.roots[1])
      # phi.limit2 = get_phase(phi.beta1.roots[2], phi.beta2.roots[2])
      
      #change to adjusted phi limits
      phi.limits = get_phaseForCI(b1.x = beta1.hat, b2.x = beta2.hat, 
                                  b1.r1 = phi.beta1.roots[1], b2.r1 = phi.beta2.roots[1], 
                                  b1.r2 = phi.beta1.roots[2], b2.r2 = phi.beta2.roots[2])
      phi.lower.limit = phi.limits$phi.lower.limit
      phi.upper.limit = phi.limits$phi.upper.limit
    }
  }else{
    phi.lower.limit = list(tan = 99, phase = 99)
    phi.upper.limit = list(tan = 99, phase = 99)
  }
  
  #calculate CI of A 
  #calculate normal ellipse parameters
  ellipse.parameters = solve.ellipse.parameters(a = B11, b = B22, c = 2*B12,
                                                d = -2*(B11*beta1.hat+B12*beta2.hat),
                                                e = -2*(B12*beta1.hat+B22*beta2.hat),
                                                f = B11*beta1.hat^2+B22*beta2.hat^2+2*B12*beta1.hat*beta2.hat-2*sigma2.hat*qf(1-alpha, 2, n-3))
  angle.point.newOrigin.major = 
    get.angle_point.newOrigin.major(x0 = ellipse.parameters$x0, 
                                    y0 = ellipse.parameters$y0, 
                                    theta.rotate = ellipse.parameters$theta.rotate)
  r.point.to.center = sqrt(ellipse.parameters$x0^2+ellipse.parameters$y0^2)
  x.new = r.point.to.center*cos(angle.point.newOrigin.major)
  y.new = r.point.to.center*sin(angle.point.newOrigin.major)
  
  #check the position of the new origin to the ellipse
  if(x.new==0&y.new==0){
    A.limit1 = ellipse.parameters$minor
    A.limit2 = ellipse.parameters$major
  }else if(x.new==0){
    A.limit1 = abs(ellipse.parameters$minor-y.new)
    A.limit2 = ellipse.parameters$minor+y.new
  }else if(y.new==0){
    A.limit1 = abs(ellipse.parameters$major-x.new)
    A.limit2 = ellipse.parameters$major+x.new
  }else if(x.new!=0&y.new!=0){
    x.new = abs(x.new)
    y.new = abs(y.new)
    fun.t = function(t){
      (ellipse.parameters$major*x.new/(t+ellipse.parameters$major^2))^2+
        (ellipse.parameters$minor*y.new/(t+ellipse.parameters$minor^2))^2-1
    }  
    
    # root.upper = 50
    # while(fun.t(-ellipse.parameters$minor^2+0.0001)*fun.t(root.upper)>0){
    #   root.upper = root.upper+50
    # }
#    root1 = uniroot(fun.t, c(-ellipse.parameters$minor^2, root.upper),extendInt="downX")$root
    root1 = uniroot(fun.t, c(-ellipse.parameters$minor^2, -ellipse.parameters$minor^2+50),extendInt="downX")$root
    x.root1 = ellipse.parameters$major^2*x.new/(root1+ellipse.parameters$major^2)
    y.root1 = ellipse.parameters$minor^2*y.new/(root1+ellipse.parameters$minor^2)
    dmin = sqrt((x.new-x.root1)^2+(y.new-y.root1)^2)
    
    # root.lower = -50
    # while(fun.t(-ellipse.parameters$major^2-0.0001)*fun.t(root.lower)>0){
    #   root.lower = root.lower-50
    # }
    #root2 = uniroot(fun.t, c(root.lower, -ellipse.parameters$major^2-0.0001), extendInt="yes")$root
    root2 = uniroot(fun.t, c(-ellipse.parameters$major^2-50, -ellipse.parameters$major^2), extendInt="upX")$root
    # root2 = uniroot(fun.t, c(-ellipse.parameters$major^2-50, -ellipse.parameters$major^2), extendInt="yes")$root
    # root2 = uniroot(fun.t, c(-ellipse.parameters$major^2-50, -ellipse.parameters$major^2), extendInt="no")$root
    x.root2 = ellipse.parameters$major^2*x.new/(root2+ellipse.parameters$major^2)
    y.root2 = ellipse.parameters$minor^2*y.new/(root2+ellipse.parameters$minor^2)
    dmax = sqrt((x.new-x.root2)^2+(y.new-y.root2)^2)
    
    A.limit1 = min(dmin, dmax)
    A.limit2 = max(dmin, dmax)
  }
  
  if(zero.in.ellipse){
    A.limit1 = 0
  }

  
  #output 
  out = list(M = list(est = m.hat, 
                      CI = c(m.hat-CI.m.hat.radius, m.hat+CI.m.hat.radius)), 
             A = list(est = A.hat, 
                      sd = se.A.hat, 
                      CI_temp = c(A.hat-qt(1-alpha/2, n-3)*se.A.hat, A.hat+qt(1-alpha/2, n-3)*se.A.hat), 
                      CI_conserve = c(A.limit1, A.limit2)), 
             phase = list(est = phase.hat, 
                          sd = se.phase.hat, 
                          CI_temp = c(phase.hat-qt(1-alpha/2, n-3)*se.phase.hat, phase.hat+qt(1-alpha/2, n-3)*se.phase.hat),
                          tan = phase.res$tan,
                          CI_tan = c(phi.lower.limit$tan, phi.upper.limit$tan), 
                          CI_conserve = c(phi.lower.limit$phase, phi.upper.limit$phase)), 
             test = list(Fstat = Fstat,
                         pval = pval,
                         R2 = MSS/TSS, 
                         acptNULL = zero.in.ellipse))
  return(out)
  
}

get_phase = function(b1.x = beta1.hat, b2.x = beta2.hat){
  ph.x = atan(-b2.x/b1.x)
  #adjust ph.x
  if(b2.x>0){
    if(ph.x<0){
      ph.x = ph.x+2*pi
    }else if (ph.x>0){
      ph.x = ph.x+pi
    }
  }else if(b2.x<0){
    if(ph.x<0){
      ph.x = ph.x+pi
    }
  }else{
    ph.x = 88#88 means one the the beta estimate is 0. I did not account for such senerio because it is rare and complicated
  }
  return(list(phase = ph.x, 
              tan = -b2.x/b1.x))
}


get_phaseForCI = function(b1.x = beta1.hat, b2.x = beta2.hat, 
                          b1.r1 = phi.beta1.roots[1], b2.r1 = phi.beta2.roots[1], 
                          b1.r2 = phi.beta1.roots[2], b2.r2 = phi.beta2.roots[2]){
  # b1.x = beta1.hat; b2.x = beta2.hat;
  # b1.r1 = phi.beta1.roots[1]; b2.r1 = phi.beta2.roots[1];
  # b1.r2 = phi.beta1.roots[2]; b2.r2 = phi.beta2.roots[2]
  
  #b1.r1 is the beta1 estimate of root 1. 
  #step 1: find which root is in the same quadrant as OLS estimate b1.x, b2.x

  x.quad = get_quad(b1.q = b1.x, b2.q = b2.x)
  phase.res = get_phase(b1.x, b2.x)
  r1.quad = get_quad(b1.q = b1.r1, b2.q = b2.r1)
  r2.quad = get_quad(b1.q = b1.r2, b2.q = b2.r2)
  
  ##the easy scenario: three are in the same quadrant 
  if(r1.quad==x.quad&r2.quad==x.quad){
    phi.limit1 = get_phase(b1.r1, b2.r1)
    phi.limit2 = get_phase(b1.r2, b2.r2)
    if(phi.limit1$phase<phi.limit2$phase){
      phi.lower.limit = phi.limit1
      phi.upper.limit = phi.limit2
    }else{
      phi.lower.limit = phi.limit2
      phi.upper.limit = phi.limit1
    }
  }else if(r1.quad == x.quad|r2.quad==x.quad){
    #more complicated scenario: one of the root is in the other quadrant
    if(r1.quad == x.quad){
      b1.s1 = b1.r1; b2.s1 = b2.r1 #s1 stands for the root that is in the same quadrant
      b1.s2 = b1.r2; b2.s2 = b2.r2
      quad.s1 = r1.quad; quad.s2 = r2.quad
    }else{
      b1.s1 = b1.r2; b2.s1 = b2.r2 #s1 stands for the root that is in the same quadrant
      b1.s2 = b1.r1; b2.s2 = b2.r1
      quad.s1 = r2.quad; quad.s2 = r1.quad
    }
    phi.s1 = get_phase(b1.s1, b2.s1)
    if(phi.s1$phase<phase.res$phase){
      phi.lower.limit = phi.s1
      s2.type = "s1 lower, s2 upper"
      phi.upper.limit = get_phaseForCI_QuadTab(s1.quad = quad.s1, s2.quad = quad.s2, s1s2 = s2.type, 
                                               b1.s2.q = b1.s2, b2.s2.q = b2.s2)
    }else{
      phi.upper.limit = phi.s1
      s2.type = "s1 upper, s1 lower"
      phi.lower.limit = get_phaseForCI_QuadTab(s1.quad = quad.s1, s2.quad = quad.s2, s1s2 = s2.type, 
                                               b1.s2.q = b1.s2, b2.s2.q = b2.s2)
    }
  }else{
    #this is the scenario that the ellipse covers three quadrants
    phi.limits = get_phaseForCI_QuadTab2(x.quad, r1.quad, r2.quad, 
                                         b1.t1 = b1.r1, b2.t1 = b2.r1, 
                                         b1.t2 = b1.r2, b2.t2 = b2.r2)
    phi.lower.limit = phi.limits$phi.lower.limit
    phi.upper.limit = phi.limits$phi.upper.limit
  }
    
  return(list(phi.lower.limit = phi.lower.limit, 
              phi.upper.limit = phi.upper.limit))
  
}

#make a table of the quadrants
get_quad = function(b1.q = b1.x, b2.q = b2.x){
  quad.table = as.data.frame(
    list(beta1.gt.0 = c(TRUE, TRUE, FALSE, FALSE), 
         beta2.gt.0 = c(TRUE, FALSE, TRUE, FALSE), 
         quad = c(1, 4, 2, 3))
  )
  quad = quad.table[(quad.table$beta1.gt.0==(b1.q>0))&(quad.table$beta2.gt.0==(b2.q>0)),
                    "quad"]
  return(quad)
}

# s1.quad = quad.s1; s2.quad = quad.s2; s1s2 = s2.type; 
# b1.s2.q = b1.s2; b2.s2.q = b2.s2
get_phaseForCI_QuadTab = function(s1.quad = quad.s1, s2.quad = quad.s2, s1s2 = s2.type, 
                                  b1.s2.q = b1.s2, b2.s2.q = b2.s2){
  # s1.quad = quad.s1; s2.quad = quad.s2; s1s2 = s2.type;
  # b1.s2.q = b1.s2; b2.s2.q = b2.s2
  quad.table2 = as.data.frame(
    list(s1.quad = c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4), 
         s2.quad = c(2, 2, 3, 3, 4, 4, 1, 1, 3, 3, 4, 4, 1, 1, 2, 2, 4, 4, 1, 1, 2, 2, 3, 3), 
         s1s2 = rep(c("s1 lower, s2 upper", "s1 upper, s1 lower"), 12), 
         s2.low = c(3, 1, 5/2, 1/2, 2, 0, 3/2, -1/2, 5/2, 1/2, 2, 0, 3/2, -1/2, 1, -1, 2, 0, 3/2, -1/2, 1, -1, 1/2, -3/2), 
         s2.up  = c(7/2, 3/2, 3, 1, 5/2, 1/2, 2, 0, 3, 1, 5/2, 1/2, 2, 0, 3/2, -1/2, 5/2, 1/2, 2, 0, 3/2, -1/2, 1, -1))
  )
  
  s2.low = quad.table2[quad.table2$s1.quad==s1.quad&quad.table2$s2.quad==s2.quad&quad.table2$s1s2==s1s2,
                       "s2.low"]
  s2.up  = quad.table2[quad.table2$s1.quad==s1.quad&quad.table2$s2.quad==s2.quad&quad.table2$s1s2==s1s2,
                       "s2.up"]
  
  phi.s2 = atan(-b2.s2.q/b1.s2.q)
  if(phi.s2<s2.low*pi){
    while(phi.s2<s2.low*pi){
      phi.s2 = phi.s2+pi
    }
  }else if(phi.s2>s2.up*pi){
    while(phi.s2>s2.up*pi){
      phi.s2 = phi.s2-pi
    }
  }
  
  if((phi.s2<(s2.up*pi))&(phi.s2>(s2.low*pi))){
    #print("#S2 is in interval")
  }else{
    print("S2 is not in the int!!!!!!!!!!!!!!!!")
  }
  return(list(phase = phi.s2, 
              tan = -b2.s2.q/b1.s2.q))
}

get_phaseForCI_QuadTab2 = function(x.quad, r1.quad, r2.quad, 
                                   b1.t1 = b1.r1, b2.t1 = b2.r1, 
                                   b1.t2 = b1.r2, b2.t2 = b2.r2){
  #t1 stands for three-quadrants, solution 1
  #This function is for when the ellipse covers three quadrants
  quad.table3 = as.data.frame(
    list(x.quad = c(1, 1, 2, 2, 3, 3, 4, 4), 
         t.quad = c(2, 4, 3, 1, 4, 2, 1, 3), 
         t.low = c(1, 2, 1/2, 3/2, 0, 1, -1/2, 1/2), 
         t.up = c(3/2, 5/2, 1, 2, 1/2, 3/2, 0, 1))
  )
  t1.low = quad.table3[quad.table3$x.quad==x.quad&quad.table3$t.quad==r1.quad, "t.low"]
  t1.up = quad.table3[quad.table3$x.quad==x.quad&quad.table3$t.quad==r1.quad, "t.up"]
  t2.low = quad.table3[quad.table3$x.quad==x.quad&quad.table3$t.quad==r2.quad, "t.low"]
  t2.up = quad.table3[quad.table3$x.quad==x.quad&quad.table3$t.quad==r2.quad, "t.up"]
  
  phi.t1 = atan(-b2.t1/b1.t1)
  if(phi.t1<t1.low*pi){
    while(phi.t1<t1.low*pi){
      phi.t1 = phi.t1+pi
    }
  }else if(phi.t1>t1.up*pi){
    while(phi.t1>t1.up*pi){
      phi.t1 = phi.t1-pi
    }
  }
  
  if((phi.t1<(t1.up*pi))&(phi.t1>(t1.low*pi))){
    #print("#t1 is in interval")
  }else{
    print("t1 is not in the int!!!!!!!!!!!!!!!!")
  }
  
  phi.t2 = atan(-b2.t2/b1.t2)
  if(phi.t2<t2.low*pi){
    while(phi.t2<t2.low*pi){
      phi.t2 = phi.t2+pi
    }
  }else if(phi.t2>t2.up*pi){
    while(phi.t2>t2.up*pi){
      phi.t2 = phi.t2-pi
    }
  }
  
  if((phi.t2<(t2.up*pi))&(phi.t2>(t2.low*pi))){
    #print("#t2 is in interval")
  }else{
    print("t2 is not in the int!!!!!!!!!!!!!!!!")
  }
  
  if(phi.t1<phi.t2){
    phi.lower.limit = list(phase = phi.t1, 
                           tan = -b2.t1/b1.t1)
    phi.upper.limit = list(phase = phi.t2, 
                           tan = -b2.t2/b1.t2)
  }else{
    phi.lower.limit = list(phase = phi.t2, 
                           tan = -b2.t2/b1.t2)
    phi.upper.limit = list(phase = phi.t1, 
                           tan = -b2.t1/b1.t1)
  }
  return(list(phi.lower.limit = phi.lower.limit, 
              phi.upper.limit = phi.upper.limit))
}

solve.ellipse.parameters = function(a = a, b = b, c = c, d = d, e = e, f = f){
  #a * x ^ 2 + b * y ^ 2 + c * x * y + d * x + e * y + f = 0
  delta1 = c^2 -4*a*b
  if(delta1>=0){
    return("Not a proper epplise")
  }else{
    #the transformation function is from wikepedia
    denominator.factor1 = 2*(a*e^2+b*d^2-c*d*e+delta1*f)
    major.minor.square.diff = sqrt(((a-b)^2+c^2))
    denominator.factor2 = c(a+b+major.minor.square.diff, a+b-major.minor.square.diff)
    major.minor = -sqrt(denominator.factor1*denominator.factor2)/delta1
    x0 = (2*b*d-c*e)/delta1
    y0 = (2*a*e-c*d)/delta1
    if(c==0){
      if(a<b){
        theta = 0
        tan.rotate = 0
      }else{
        theta = pi/2
        tan.rotate = Inf
      }
    }else{
      tan.rotate = (b-a-major.minor.square.diff)/c
      theta = atan(tan.rotate)
    }
    return(list(major = major.minor[1], 
                minor = major.minor[2], 
                x0 = x0, 
                y0 = y0, 
                theta.rotate = theta, 
                tan.rotate = tan.rotate))
  }
}

get.angle_point.newOrigin.major = function(x0 = ellipse.parameters$x0, 
                                           y0 = ellipse.parameters$y0, 
                                           theta.rotate = ellipse.parameters$theta.rotate){
  #if both x0=0 and y0=0, then the distance is already there:
  #min distance = b; max distance = a
  if(x0==0&y0==0){
    return("ellipse is centered at (0, 0)")
  }else if(y0==0){
    angle = theta.rotate
    return(angle)
  }else if(x0==0){
    angle = pi/2-theta.rotate
    return(angle)
  }else{
    tan.center.0.x_axis = y0/x0
    theta.center.0.x_axis = atan(tan.center.0.x_axis)
    if(theta.center.0.x_axis*theta.rotate<0){
      angle = abs(theta.center.0.x_axis)+abs(theta.rotate)
      return(angle)
    }else if(theta.center.0.x_axis*theta.rotate>=0){
      angle = abs(theta.center.0.x_axis-theta.rotate)
      return(angle)
    }
  }
}


##END OF COSINOR CODE--------------



#SET UP DRAWING CIRCADIAN PLOTS
circadianDrawing <- function(tod, expr, apar, labels, specInfo=NULL){	
  getPred <- function(parS, xx) {	
    parS$A * cos(2*pi/24 * xx + parS$phase) + parS$offset
  }
  peak <- round(apar$peak)
  amain <- paste0(' peak = ',peak, " (pval=", round(apar$pvalue, 4), ")",
                  "R2 = ", round(apar$R2, 2))
  times <- seq(-7,18,0.1)
  pred <- getPred(apar,times)
  
  labelColor <- as.numeric(factor(labels))
  
  plot(tod,expr,col=labelColor, pch=16,cex=2,
       main=amain,xlim=c(-7,18),
       xaxt="n", yaxt="n",
       xlab='',ylab='')
  ytick = pretty(par("usr")[3:4])
  yl = formatC(ytick, format="f", digits=2) 
  axis(2,cex.axis=2.2,at=ytick,labels = yl)
  mtext("Expression", side=2, line=2.6, cex=2.2)
  xtick = pretty(par("usr")[1:2])
  xl = formatC(xtick, format="f", digits=0) 
  axis(1,cex.axis=2.2,at=xtick,labels = xl)
  mtext("TOD", side=1, line=2.6, cex=2.2)
  
  smoothingSpline = smooth.spline(times, pred, spar=0.35)
  lines(smoothingSpline,col='red',lwd=4)
  
  box(which = "plot", lty = "solid",lwd=3)	
}

#END CIRCADIAN PLOTS



#ENTER DATA AND SET UP FOR RUNNING ANALYSIS

#logcpm values format: column names are subject ids and row names are gene ids
data <- read.csv("~/Desktop/Panos Collab/v0_pseudoExpr/v0_Astro.csv", row.names = 1)
#TOD values and matching subject IDs
pheno <- read.csv ("~/Desktop/Panos Collab/pb_metadata.csv", row.names = 1)

#For the single cell pilot I did a bunch of filtering and subsetting - probably will not be necessary here but just in case:

#add in columns with "1" for each cell type - filter by cell type 
pheno_cell <- filter(pheno, pheno$Astro == "1")

#filter by group - NA = Non-AD, AD = Alzheimer's, Com = Small comparison group (no SZ)
pheno_cell_NA <- filter(pheno_cell, pheno_cell$Dx == "Control")
pheno_cell_AD <- filter(pheno_cell, pheno_cell$Dx == "AD")
pheno_cell_Com <- filter(pheno_cell, pheno_cell$COM == "1")

#THIS IS PROBABLY THE PART WE NEED 
#Make a data table- sub.data = data[,match(pheno$ID,colnames(data))]
cell.data = data[,match(pheno_cell$ID,colnames(data))]
colnames(data)

#Filter by group 
cell_Com.data = cell.data[,which((colnames(cell.data) %in% pheno_cell_Com$ID) == TRUE)]

#Way to check your names - all(pheno_Astro$ID == colnames(data))
all(pheno_cell_Com$ID == colnames(cell_Com.data))
SubID <- pheno_Astro_Com$SubID

#END DATA SET UP 


#ACTUAL RUNNING OF THE COSINOR CODE
out.list = mclapply(1:nrow(cell_Com.data), function(i){ 
  if(i%%1000==0) print(i)
  res = one_cosinor_OLS(tod = as.numeric(pheno_cell_Com$TOD), y = as.numeric(cell_Com.data[i,]))
  res.onerow = as.data.frame(list(offset = res$M$est, offset.ll = res$M$CI[1], offset.ul = res$M$CI[2], 
                                  A = res$A$est, A.sd = res$A$sd, 
                                  phase = res$phase$est, phase.sd = res$phase$sd, 
                                  pvalue = res$test$pval, R2 = res$test$R2))
  return(res.onerow)
}, mc.cores = 1) 

out.tab = do.call(rbind.data.frame, out.list)
out.tab$qvalue = p.adjust(out.tab$pvalue, "BH")
out.tab$peak = 24 - out.tab$phase*24/(2*pi)
row.names(out.tab) = row.names(cell.data)
setwd("~/Desktop/Panos Collab/v0_rhythmicity/Astro//")
write.csv(out.tab,paste0("observed_Astro_Com.csv"))



#CIRCADIAN PLOTS OF CORE GENES - IF ERROR TRY INDIVIDUALLY 
tod = pheno_cell_AD$TOD
#specInfo = paste0(sex)
labels <-rep(1, ncol(cell_AD.data))

core_symbols = c("ARNTL", "CLOCK", "NPAS2", "PER1", "PER2", "PER3", "CRY1", "CRY2", "RORA", "RORB", "NR1D1", "NR1D2", "CIART", "OPRL1", "DBP")

agene <- core_symbols[i]
for(i in 1:length(agene)){
  fileName <- paste0("_",agene,'.pdf')
  pdf(fileName)
  circadianDrawing(tod=tod, expr=unlist(cell_AD.data[agene,]), apar=out.tab[agene,],labels=labels)
  dev.off()
}

