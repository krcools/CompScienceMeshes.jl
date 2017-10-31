export gauss_quad_linear,gauss_points
function gauss_quad_linear(N,a,b)
    N=N-1;
    N1=N+1;
    N2=N+2;
    xu=linspace(-1,1,N1)';
    # Initial guess
    y=cos.((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin.(pi*xu*N/N2);
    # Legendre-Gauss Vandermonde Matrix
    L=zeros(N1,N2);
    # Derivative of LGVM
    Lp=zeros(N1,N2);
    # Compute the zeros of the N+1 Legendre Polynomial
    # using the recursion relation and the Newton-Raphson method
    y0=2;
    # Iterate until new points are uniformly within epsilon of old points
    while (maximum(abs.(y-y0))>eps(1.0))
        L[:,1]=1;
        Lp[:,1]=0;
        L[:,2]=y;
        #Lp[:,2]=1;
        for k=2:N1
            L[:,k+1]=( (2*k-1)*y'.*L[:,k]-(k-1)*L[:,k-1] )/k;
        end
        Lp=(N2)*( L[:,N1]-y'.*L[:,N2] )./(1-y.^2)';
        y0=y;
        y=y0'-L[:,N2]./Lp;
        y=y'
    end
    y=y'
    # Linear map from[-1,1] to [a,b]
    x=(a*(1-y)+b*(1+y))/2;
    # Compute the weights
    w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;
    return x, w
end

function gauss_points(N)
    x,w = gauss_quad_linear(N,-1,1)
    gpt = zeros(N*N,2)
    gpw = zeros(N*N,1)
    for i = 1:N
        u = x[i]
        wu = w[i]
        for j = 1:N
            v = x[j]
            wv = w[j]
            gpt[(i-1)*N + j,1] = u
            gpt[(i-1)*N + j,2] = v
            gpw[(i-1)*N + j] = wu * wv
        end
    end
    return gpt,gpw
end
