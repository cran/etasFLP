\name{etasFLP-internal}
\title{Internal etasFLP functions}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	elementary functions to solve a triangle given three elements
%	a,b,c sides; alpha,beta gamma angles (in radians) 
%	alpha opposite to a; beta opposite to b; gamma opposite to c 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\alias{triangle.alphabetac.b}
\alias{triangle.abgamma.c}
\alias{triangle.abc.gamma}
\alias{triangle.abgamma.alpha}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	embedding.rect.cat: returns the rectangular region 
%	embedding latitude and lomngitude of the points (x,y) 
%	of a catalog

%	embedding.rect.cat.eps returns an enlarged region 
%	embedding latitude and lomngitude of the points (x,y) 
%	of a catalogs to avoid to have points on the boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\alias{embedding.rect.eps}
\alias{embedding.rect.cat.eps}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% other undocumented functions

% angle.positive	force an angle in -2pi-2pi to 0-2pi

\alias{angle.positive}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	cartesian2polar2D: from cartesian to polar coordinates (in 2 D) 
%
%	x is a matrix with n rows and two columns with the cartesian coordinates 
%        of n points (column1=x, column2=y)
%
%	orig is the point of the new origin
%	
%	returned value a matrix of n rows and two columns with polar coordinates
%	(column1=rho, column2=theta)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\alias{cartesian2polar2D}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  region.P
%%%	the region should be defined NOT cicly (n different points)
%%%	transform to polar coordinates region from an origin P (inside the region)
%%%	and then compute the length of k  segments from P to the border of the region
%%%  (the segments start from P at angles of 2pi/k)
%%%
%%% written in 2006-2007
%% comments added in 2011
%%
%%
%% region can be a generic polygonal region (convex only!)
%%
%% polar1_graphichelp.pdf contains an help picture
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\alias{region.P}


\alias{flpkspace}
\alias{etas.mod2}
% \alias{etasflp.integral} dropped beginning from version 1.3.0
\alias{flp1.etas.nlm}
\alias{longlat.to.km}
\alias{timechar.to.seq}

\description{Internal functions used in the package etasFLP }
\author{Marcello Chiodi}
\details{These are not intended to be used directly by users}
\keyword{internal}