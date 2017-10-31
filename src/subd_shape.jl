export shape_function,shape_function_der,shape_function_der2
function boxspline_func(v,w)
    u=1.0-v-w;

    x = zeros(1,12)

    x[1] = (u*u*u*u + 2.0*u*u*u*v)/12.0

    x[2] = (u*u*u*u + 2.0*u*u*u*w)/12.0

    x[3] = (u*u*u*u + 2.0*u*u*u*w + 6.0*u*u*u*v + 6.0*u*u*v*w + 12.0*u*u*v*v + 6.0*u*v*v*w + 6.0*u*v*v*v + 2.0*v*v*v*w + v*v*v*v) / 12.0

    x[4] = (6.0*u*u*u*u + 24.0*u*u*u*w + 24.0*u*u*w*w + 8.0*u*w*w*w + w*w*w*w + 24.0*u*u*u*v + 60.0*u*u*v*w + 36.0*u*v*w*w + 6.0*v*w*w*w + 24.0*u*u*v*v + 36.0*u*v*v*w + 12.0*v*v*w*w + 8.0*u*v*v*v + 6.0*v*v*v*w + v*v*v*v) / 12.0

    x[5] = (u*u*u*u + 6.0*u*u*u*w + 12.0*u*u*w*w + 6.0*u*w*w*w + w*w*w*w + 2.0*u*u*u*v + 6.0*u*u*v*w + 6.0*u*v*w*w + 2.0*v*w*w*w) / 12.0

    x[6] = (2.0*u*v*v*v + v*v*v*v) / 12.0

    x[7] = (u*u*u*u + 6.0*u*u*u*w + 12.0*u*u*w*w + 6.0*u*w*w*w + w*w*w*w + 8.0*u*u*u*v + 36.0*u*u*v*w + 36.0*u*v*w*w + 8.0*v*w*w*w + 24.0*u*u*v*v + 60.0*u*v*v*w + 24.0*v*v*w*w + 24.0*u*v*v*v + 24.0*v*v*v*w + 6.0*v*v*v*v) / 12.0

    x[8] = (u*u*u*u + 8.0*u*u*u*w + 24.0*u*u*w*w + 24.0*u*w*w*w + 6.0*w*w*w*w + 6.0*u*u*u*v + 36.0*u*u*v*w + 60.0*u*v*w*w + 24.0*v*w*w*w + 12.0*u*u*v*v + 36.0*u*v*v*w + 24.0*v*v*w*w + 6.0*u*v*v*v + 8.0*v*v*v*w + v*v*v*v)/12.0

    x[9] = (2.0*u*w*w*w + w*w*w*w) / 12.0

    x[10] = (2.0*v*v*v*w + v*v*v*v) / 12.0

    x[11] = (2.0*u*w*w*w + w*w*w*w + 6.0*u*v*w*w + 6.0*v*w*w*w + 6.0*u*v*v*w + 12.0*v*v*w*w + 2.0*u*v*v*v + 6.0*v*v*v*w + v*v*v*v) / 12.0

    x[12] = (w*w*w*w + 2.0*v*w*w*w) / 12.0

    return x'
end

function boxspline_func_der(v,w)
    u=1.0-v-w
    der1=zeros(2,12)
    der1[1,1] = (-6.0*v*pow(u,2.0) - 2.0*pow(u,3.0))/12.0
    der1[2,1] = (-6.0*v*pow(u,2.0) - 4.0*pow(u,3.0))/12.0

    der1[1,2] = (-4.0*pow(u,3.0)-6.0*pow(u,2.0)*w)/12.0
    der1[2,2] = (-2.0*pow(u,3.0)-6.0*pow(u,2.0)*w)/12.0

    der1[1,3] = (-2.0*pow(v,3.0)-6.0*pow(v,2.0)*u + 6.0*v*pow(u,2.0)+2.0*pow(u,3.0))/12.0
    der1[2,3] = (-4.0*pow(v,3.0)-18.0*pow(v,2.0)*u - 12.0*v*pow(u,2.0)-2.0*pow(u,3.0) - 6.0*pow(v,2.0)*w-12.0*v*u*w - 6.0*pow(u,2.0)*w)/12.0

    der1[1,4] = (-4.0*pow(v,3.0)-24.0*pow(v,2.0)*u - 24.0*v*pow(u,2.0)-18.0*pow(v,2.0)*w - 48.0*v*u*w-12.0*pow(u,2.0)*w - 12.0*v*pow(w,2.0) - 12.0*u*pow(w,2.0) - 2.0*pow(w,3.0))/12.0
    der1[2,4] = (-2.0*pow(v,3.0)-12.0*pow(v,2.0)*u - 12.0*v*pow(u,2.0)-12.0*pow(v,2.0)*w - 48.0*v*u*w-24.0*pow(u,2.0)*w - 18.0*v*pow(w,2.0)-24.0*u*pow(w,2.0) - 4.0*pow(w,3.0))/12.0

    der1[1,5] = (-6.0*v*pow(u,2.0)-2.0*pow(u,3.0) - 12.0*v*u*w-12.0*pow(u,2.0)*w - 6.0*v*pow(w,2.0)-18.0*u*pow(w,2.0) - 4.0*pow(w,3.0))/12.0

    der1[2,5] = (2.0*pow(u,3.0)+6.0*pow(u,2.0)*w - 6.0*u*pow(w,2.0)-2.0*pow(w,3.0))/12.0

    der1[1,6] = (2.0*pow(v,3.0)+6.0*pow(v,2.0)*u)/12.0
    der1[2,6] = -pow(v,3.0)/6.0

    der1[1,7] = (24.0*pow(v,2.0)*u+24.0*v*pow(u,2.0) + 4.0*pow(u,3.0)+12.0*pow(v,2.0)*w + 48.0*v*u*w+18.0*pow(u,2.0)*w + 12.0*v*pow(w,2.0)+12.0*u*pow(w,2.0) + 2.0*pow(w,3.0))/12.0
    der1[2,7] = (12.0*pow(v,2.0)*u+12.0*v*pow(u,2.0) + 2.0*pow(u,3.0)-12.0*pow(v,2.0)*w + 6.0*pow(u,2.0)*w-12.0*v*pow(w,2.0) - 6.0*u*pow(w,2.0)-2.0*pow(w,3.0))/12.0

    der1[1,8] = (-2.0*pow(v,3.0)-6.0*pow(v,2.0)*u + 6.0*v*pow(u,2.0)+2.0*pow(u,3.0) - 12.0*pow(v,2.0)*w+12.0*pow(u,2.0)*w - 12.0*v*pow(w,2.0)+12.0*(1.0-v-w)*pow(w,2.0))/12.0
    der1[2,8] = (2.0*pow(v,3.0)+12.0*pow(v,2.0)*u + 18.0*v*pow(u,2.0)+4.0*pow(u,3.0) + 12.0*pow(v,2.0)*w+48.0*v*u*w + 24.0*pow(u,2.0)*w+12.0*v*pow(w,2.0) + 24.0*u*pow(w,2.0))/12.0

    der1[1,9] = -pow(w,3.0)/6.0
    der1[2,9] = (6.0*u*pow(w,2.0)+2.0*pow(w,3.0))/12.0

    der1[1,10] = (4.0*pow(v,3.0)+6.0*pow(v,2.0)*w)/12.0
    der1[2,10] = pow(v,3.0)/6.0

    der1[1,11]= (2.0*pow(v,3.0)+6.0*pow(v,2.0)*u + 12.0*pow(v,2.0)*w+12.0*v*u*w + 18.0*v*pow(w,2.0)+6.0*u*pow(w,2.0) + 4.0*pow(w,3.0))/12.0

    der1[2,11]= (4.0*pow(v,3.0)+6.0*pow(v,2.0)*u + 18.0*pow(v,2.0)*w+12.0*v*u*w + 12.0*v*pow(w,2.0)+6.0*u*pow(w,2.0) + 2.0*pow(w,3.0))/12.0

    der1[1,12] = pow(w,3.0)/6.0
    der1[2,12] = (6.0*v*pow(w,2.0)+4.0*pow(w,3.0))/12.0

    return der1'

end

function boxspline_func_der2(v,w)
    u = 1.0-v-w
    der2=zeros(3,12)
    der2[1,1] = v*u
    der2[2,1] = v*u+pow(u,2.0)
    der2[3,1] = (12.0*v*u+6.0*pow(u,2.0))/12.0

    der2[1,2] = pow(u,2.0)+u*w
    der2[2,2] = u*w
    der2[3,2] = (6.0*pow(u,2.0)+12.0*u*w)/12.0

    der2[1,3] = -2.0*v*u
    der2[2,3] = pow(v,2.0)+v*u+v*w+u*w
    der2[3,3] = (6.0*pow(v,2.0)-12.0*v*u -6.0*pow(u,2.0))/12.0

    der2[1,4] = pow(v,2.0)-2.0*pow(u,2.0) + v*w-2.0*u*w
    der2[2,4] = -2.0*v*u-2.0*pow(u,2.0) + v*w+pow(w,2.0)
    der2[3,4] = (6.0*pow(v,2.0)-12.0*pow(u,2.0) + 24.0*v*w+6.0*pow(w,2.0))/12.0

    der2[1,5] = v*u + v*w + u*w + pow(w,2.0)
    der2[2,5] = - 2.0*u*w
    der2[3,5] = (- 6.0*pow(u,2.0) - 12.0*u*w + 6.0*pow(w,2.0))/12.0

    der2[1,6] = v*u
    der2[2,6] = 0.0
    der2[3,6] = -pow(v,2.0)/2.0

    der2[1,7] = (-24.0*pow(v,2.0)+12.0*pow(u,2.0)-24.0*v*w + 12.0*u*w)/12.0
    der2[2,7] = (-24.0*pow(v,2.0)-24.0*v*u-24.0*v*w - 24.0*u*w)/12.0
    der2[3,7] = (-12.0*pow(v,2.0)+6.0*pow(u,2.0)-24.0*v*w - 12.0*u*w-6.0*pow(w,2.0))/12.0

    der2[1,8] = -2.0*v*u-2.0*v*w-2.0*u*w - 2.0*pow(w,2.0)
    der2[2,8] = v*u+pow(u,2.0)-2.0*v*w - 2.0*pow(w,2.0)
    der2[3,8] = (-6.0*pow(v,2.0)-12.0*v*u+6.0*pow(u,2.0) - 24.0*v*w-12.0*pow(w,2.0))/12.0

    der2[1,9] = 0.0
    der2[2,9] = u*w
    der2[3,9] = -pow(w,2.0)/2.0

    der2[1,10] = (12.0*pow(v,2.0)+12.0*v*w)/12.0
    der2[2,10] = 0.0
    der2[3,10] = pow(v,2.0)/2.0

    der2[1,11]= (12.0*v*u+12.0*v*w+12.0*u*w + 12.0*pow(w,2.0))/12.0
    der2[2,11]= pow(v,2.0)+v*u+v*w+u*w
    der2[3,11]= (6.0*pow(v,2.0)+12.0*v*u+24.0*v*w + 12.0*u*w+6.0*pow(w,2.0))/12.0

    der2[1,12]= 0.0
    der2[2,12]= v*w+pow(w,2.0)
    der2[3,12]= pow(w,2.0)/2.0

  return der2';

end

function reg_shapefun(v,w,d)
    m = pickmatrx_reg()
    if d == 0
        a = boxspline_func(v,w)
        func = m'*a
        return func
    elseif d == 1
        b = boxspline_func_der(v,w)
        der = m'*b
        return der
    elseif d == 2
        c = boxspline_func_der2(v,w)
        der3 = m'*c
        return der3
    else
        error("currently we can only obtain upto second derivative")
    end
end

function pickmatrx_reg()
    m = zeros(12,12)
    m[1,4] = 1
    m[2,5] = 1
    m[3,3] = 1
    m[4,1] = 1
    m[5,6] = 1
    m[6,10] = 1
    m[7,2] = 1
    m[8,7] = 1
    m[9,12] = 1
    m[10,9] = 1
    m[11,8] = 1
    m[12,11] = 1
    return m
end

function irreg_shapefun(v,w,ival,d)
    eps = 1.0e-10
    u = 1.0 - v - w
    na = Int64(0)
    min = 0.0
    max = 0.5
    # evaluate the number of the required subdivisions
    while !((u>(min-eps))&&(u<(max+eps)))
        na += 1
        min = max
        max += 1.0/pow(2.0,(na+1))
    end
    potenz = na+1
    pow2 = pow(2.0, na)
    v *= pow2
    w *= pow2
    u = 1.0-v-w
    # coordinate transformation
    if (v > (0.5-eps))
        v = 2.0*v-1.0
        w = 2.0*w
        pm = pickmtrx_irreg(ival,1)
        jfac = pow(2.0, na+1)
    elseif (w >(0.5-eps))
        v = 2.0*v
        w = 2.0*w-1.0
        pm = pickmtrx_irreg(ival,3)
        jfac = pow(2.0, na+1)
    else
        v = 1.0 - 2.0*v
        w = 1.0-2.0*w
        pm = pickmtrx_irreg(ival,2)
        jfac = pow(2.0, na+1)
    end
    u = 1.0-v-w
    # compute matrix A
    A,A_hat = matrix_A(ival)
    # compute p*power(A,n)
    for ip = 1:potenz-1
        A_hat*=A
    end
    An = A_hat
    PAA = An'*pm'

    if d == 0
        shape = boxspline_func(v,w)
        shape_irreg = PAA*shape
        return shape_irreg
    elseif d == 1
        shape_der = boxspline_func_der(v,w)
        shape_irreg_der = jfac * PAA * shape_der
        return shape_irreg_der
    elseif d == 2
        shape_der2 = boxspline_func_der2(v,w)
        shape_irreg_der2 = jfac * jfac * PAA * shape_der2
        return shape_irreg_der2
    else
        error("currently we can only obtain upto second derivative")
    end
end

function pow(x,y)
    return x^y
end

function pickmtrx_irreg(ival,case)
    p=zeros(12,ival+12)
    if case == 1
        p[1,3]= 1.0
        p[2,1]= 1.0
        p[3,ival+4]= 1.0
        p[4,2]= 1.0
        p[5,ival+1]  = 1.0
        p[6,ival+9]= 1.0
        p[7,ival+3]= 1.0
        p[8,ival+2]= 1.0
        p[9,ival+5]= 1.0
        p[10,ival+8]= 1.0
        p[11,ival+7]= 1.0
        p[12,ival+10]= 1.0
    elseif case == 2
        p[1,ival+10]= 1.0
        p[2,ival+7]= 1.0
        p[3,ival+5]= 1.0
        p[4,ival+2]= 1.0
        p[5,ival+3]= 1.0
        p[6,ival+6]= 1.0
        p[7,ival+1]  = 1.0
        p[8,2]     = 1.0
        p[9,ival+4]= 1.0
        p[10,ival]= 1.0
        p[11,1]    = 1.0
        p[12,3]    = 1.0
    elseif case == 3
        p[1,1]= 1.0;
        p[2,ival]  = 1.0;
        p[3,2]= 1.0;
        p[4,ival+1]    = 1.0;
        p[5,ival+6]= 1.0;
        p[6,ival+3]  = 1.0;
        p[7,ival+2]= 1.0;
        p[8,ival+5]  = 1.0;
        p[9,ival+12]= 1.0;
        p[10,ival+7]  = 1.0;
        p[11,ival+10]= 1.0;
        p[12,ival+11]= 1.0;
    else
        error("Case number wrong")
    end
    return p
end

function matrix_A(ival)
    if ival == 3
        alpha = 0.5625
        b = 3.0/16.0
    else
        alpha = 0.375
        b = 3.0/(8.0 * ival)
    end
    c = 3.0/8.0
    d = 1.0/8.0
    S = zeros(ival+1,ival+1)
    S[1,1]=1-alpha
    for i = 2:ival+1
        S[1,i] = b
    end
    S[2,1]=c
    S[2,2]=c
    S[2,3]=d
    S[2,ival+1]=d
    S[ival+1,1]=c
    S[ival+1,2]=d
    S[ival+1,ival]=d
    S[ival+1,ival+1]=c
    for i = 3:ival
        S[i,1] = c
        S[i,i-1] = d
        S[i,i] = c
        S[i,i+1] = d
    end
    s22 = [0.375  0.125   0.0   0.0   0.0;
    0.125  0.375   0.125  0.0   0.0;
    0.0    0.125   0.375  0.0   0.0;
    0.375  0.0     0.0    0.125 0.0;
    0.125  0.0     0.0    0.375 0.125;
    0.0    0.0     0.0    0.125 0.375]

    s12 = [0.125  0.0    0.0    0.0    0.0;
    0.0625 0.0625 0.0625 0.0    0.0;
    0.0    0.0    0.125  0.0    0.0;
    0.0625 0.0    0.0    0.0625 0.0625;
    0.0    0.0    0.0    0.0    0.125]

    s11 = zeros(5,ival+1)
    s21 = zeros(6,ival+1)

    s11[1,1] = 0.125
    s11[1,2] = 0.375
    s11[1,ival+1] = 0.375
    s11[2,1] = 0.0625
    s11[2,2] = 0.625
    s11[2,3] = 0.0625
    s11[2,ival+1] = 0.0625
    s11[3,1] = 0.125
    s11[3,2] = 0.375
    s11[3,3] = 0.375
    s11[4,1] = 0.0625
    s11[4,2] = 0.0625
    s11[4,ival] = 0.0625
    s11[4,ival+1] = 0.625
    s11[5,1] = 0.125
    s11[5,ival] = 0.375
    s11[5,ival+1] = 0.375

    s21[1,2] = 0.375
    s21[1,ival+1] = 0.125
    s21[2,2] = 0.375
    s21[3,2] = 0.375
    s21[3,3]  = 0.125
    s21[4,2] = 0.125
    s21[4,ival+1] = 0.375
    s21[5,ival+1] = 0.375
    s21[6,ival]= 0.125
    s21[6,ival+1]= 0.375

    A = zeros(ival+6,ival+6)
    A[1:ival+1,1:ival+1] = S
    A[ival+2:ival+6,1:ival+1] = s11
    A[ival+2:ival+6,ival+2:ival+6] = s12
    A_hat = zeros(ival+12,ival+6)
    A_hat[1:ival+6,1:ival+6] = A
    A_hat[ival+7:ival+12,1:ival+1] = s21
    A_hat[ival+7:ival+12,ival+2:ival+6] = s22

    return A,A_hat
end

function shape_function(v,w,ival)
    if ival == 6
        return reg_shapefun(v,w,0)
    else
        return irreg_shapefun(v,w,ival,0)
    end
end

function shape_function_der(v,w,ival)
    if ival == 6
        return reg_shapefun(v,w,1)
    else
        return irreg_shapefun(v,w,ival,1)
    end
end

function shape_function_der2(v,w,ival)
    if ival == 6
        return reg_shapefun(v,w,2)
    else
        return irreg_shapefun(v,w,ival,2)
    end
end
