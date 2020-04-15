def gompertz(t, a, b, c):
    
    ct = c*t
    
    v = np.exp(-ct)
    g = b*v
    
    u = np.exp(-g)
    f = a*u
    
    r = c*g
    
    w = f*v
    
    
    da = u
    db = -w
    dc = f*g*t
    
    d11 = np.zeros_like(t)
    d12 = db/a
    d13 = dc/a
    
    d22 = w*v
    d23 = w*t*(1.-g)
    d33 = -d23*t*b
    
    
    df = np.array([da, db, dc])
    
    d2f = np.array([[d11, d12, d13], [d12, d22, d23], [d13, d23, d33]])
    
    
    return f, r, df, d2f
