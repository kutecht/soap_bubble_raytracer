//
//  RayTraceImageView.m
//  SoapBubble
//
//  Created by Kevin Utecht on 6/21/13.
//
//  Modeling Interference Color
//

#import "RayTracer.h"

#define MAX_RECURSION_LEVEL   3

typedef struct
{
    double x,y,z;
} vector_s;


typedef struct
{
    double r,g,b;
} color_s;


typedef struct
{
    vector_s center;
    double radius;
    double yMin;
    double yMax;
    double interferencePercentage;
    double reflectionPercentage;
    color_s reflection;
    color_s ambient;
    color_s diffuse;
    color_s specular;
    double specularWidthFactor;
    double reflectionDiffusion;
    vector_s timeMotionCoefficient;
    double thicknessCoefficient;
} bubble_s;


typedef struct
{
    vector_s orgin;
    double angle;
} light_s;


typedef struct
{
    vector_s origin;
    vector_s direction;
    bubble_s *source;
    double diffusionAngle;
} ray_s;


typedef struct
{
    bubble_s *bubble;
    double position;
    vector_s normalAtPoint;
} intersect_s;


// the CIE tristimulus variables
// Source: Color Science: Concepts and Methods, Quantitative Data and Formulae, second edition, Wiley, New York, 1982
float C[40] = {33.0, 47.40, 63.30, 80.60, 98.10, 112.40, 121.50, 124.00, 123.10, 123.80, 123.90, 120.70, 112.10,
    102.30, 96.90, 98.00, 102.10, 105.2, 105.30, 102.30, 97.80, 93.20, 89.70, 88.40, 88.10, 88.00,
    87.80, 88.20, 87.90, 86.30, 84.00, 80.20, 76.30, 72.40, 68.30, 64.40, 61.50, 59.20, 58.10, 58.20};

float x_bar[40] = {0.0014, 0.0042, 0.0143, 0.0435, 0.1344, 0.2839, 0.3483, 0.3362, 0.2908, 0.1954, 0.0956, 0.0320,
    0.0049, 0.0093, 0.0633, 0.1655, 0.2904, 0.4334, 0.5945, 0.7621, 0.9163, 1.0263, 1.0622, 1.0026,
    0.8544, 0.6424, 0.4479, 0.2835, 0.1649, 0.0874, 0.0468, 0.0227, 0.0114, 0.0058, 0.0029, 0.0014,
    0.0007, 0.0003, 0.0002, 0.0001};

float y_bar[40] = {0.0000, 0.0001, 0.0004, 0.0012, 0.0040, 0.0116, 0.0230, 0.0380, 0.0600, 0.0910, 0.1390, 0.2080,
    0.3230, 0.5030, 0.7100, 0.8620, 0.9540, 0.9950, 0.9950, 0.9520, 0.8700, 0.7570, 0.6310, 0.5030,
    0.3810, 0.2650, 0.1750, 0.1070, 0.0610, 0.0320, 0.0170, 0.0082, 0.0041, 0.0021, 0.0010, 0.0005,
    0.0002, 0.0001, 0.0001, 0.0000};

float z_bar[40] = {0.0065, 0.0201, 0.0679, 0.2074, 0.6456, 1.3856, 1.7471, 1.7721, 1.6692, 1.2876, 0.8130, 0.4652,
    0.2720, 0.1582, 0.0782, 0.0422, 0.0203, 0.0087, 0.0039, 0.0021, 0.0017, 0.0011, 0.0008, 0.0003,
    0.0002, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
    0.0000, 0.0000, 0.0000, 0.0000};


vector_s vadd(vector_s a, vector_s b)
{
    vector_s c;
    
    c.x = a.x + b.x;
    c.y = a.y + b.y;
    c.z = a.z + b.z;
    return c;
}


vector_s vsub(vector_s a, vector_s b)
{
    vector_s c;
    
    c.x = a.x - b.x;
    c.y = a.y - b.y;
    c.z = a.z - b.z;
    return c;
}


vector_s vneg(vector_s a)
{
    vector_s b;
    
    b.x = -a.x;
    b.y = -a.y;
    b.z = -a.z;
    return b;
}


vector_s svproduct(double k, vector_s a)
{
    a.x *= k;
    a.y *= k;
    a.z *= k;
    return a;
}


double vdot(vector_s a, vector_s b)
{
    return (a.x * b.x + a.y * b.y + a.z * b.z);
}


vector_s vcross(vector_s a, vector_s b)
{
    vector_s c;
    
    c.x = a.y * b.z - b.y * a.z;
    c.y = b.x * a.z - a.x * b.z;
    c.z = a.x * b.y - b.x * a.y;
    return c;
}


vector_s normalize(vector_s a)
{
    double len;
    vector_s b;
    
    len = sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
    b.x = a.x / len;
    b.y = a.y / len;
    b.z = a.z / len;
    return b;
}


double sqr (double x)
{
    x = x * x;
    return(x);
}


intersect_s sphere(bubble_s *bubble, ray_s r, double currentTime)
{
    vector_s v;
    vector_s normal;
    vector_s center;
    intersect_s i;
    double b, c, d;
    double sol1, sol2;
    
    i.bubble = NULL;
    
    // find out what the center is at this time
    center = svproduct(currentTime, bubble->timeMotionCoefficient);
    center = vadd(bubble->center, center);
    
    v = vsub(r.origin, center);
    b = 2 * vdot(r.direction, v);
    c = vdot(v, v) - bubble->radius;
    
    d = b * b - 4 * c;
    if( d < 0 )
        return i;
    d = sqrt(d);
    sol1 = ( -b + d ) / 2;
    sol2 = ( -b - d ) / 2;
    if( sol1 <= 0 )
        sol1 = sol2;
    if ( sol2 <= 0 )
        sol2 = sol1;
    i.position = (sol1 < sol2) ? sol1 : sol2 ;
    if (i.position <= 0)
    {
        // intersection is behind eye
        return i;
    }
    i.bubble = bubble;
    
    // normal = the direction of the radius
    normal = vsub(vadd(r.origin, svproduct(i.position, r.direction)), center);
    i.normalAtPoint = normalize(normal);
    
    return i;
}


intersect_s intersect(ray_s r, bubble_s *bubble, double currentTime)
{
    intersect_s inter, intermin;
    
    intermin.bubble = NULL;
    
    inter = sphere(bubble, r, currentTime);
    
    // update the minimum intersection distance if the new intersection exists (ray intersects the object), and the intersection distance
    // is smaller that the one logged and also the object intersected is not the object that the ray is originating from.
    if (inter.bubble && (!intermin.bubble || (inter.position < intermin.position)) && (inter.bubble != r.source))
    {
        intermin = inter;
    }
    
    return intermin;
}


double quickcos(double x)
{
    double val;
    
    val = 1 - x*x/2.4684;
    return val;
}


double quickinvcos(double x)
{
    double val;
    
    val = sqrt(2.4684*(1-x));
    return val;
}


// floating point random number generator
double rnd()
{
    double t;
    
    t = (double)(random() & 0xffffff) / 16777215;
    return t;
}


// random number between -1 and 1
double rand1()
{
    return ((rnd() - .5) * 2);
}


// approximate a gaussian random number generator between -1 and 1
double grand()
{
    double t;
    
    t = rand1();
    
    // and then square it.  Don't lose the sign!
    return (t*ABS(t));
}


// pick a random ray somewhere inside the solid angle
ray_s sample_ray(ray_s r, double theta)
{
    double phi1, phi2;
    vector_s c; // directional cosines
    ray_s r2;
    
    // adds the following 2 angles in 2 of the 3 angles specified by the
    // directional cosines
    phi1 = grand() * theta;
    phi2 = grand() * theta;
    
    c = r.direction;
    
    // choose two of them that are linearly independent
    if (c.x > .5)
    {
        c.x = quickcos(quickinvcos(c.x) + phi1);
        c.y = quickcos(quickinvcos(c.y) + phi2);
    }
    else
    {
        if (c.y < .5)
        {
            c.y = quickcos(quickinvcos(c.y) + phi1);
            c.z = quickcos(quickinvcos(c.z) + phi2);
        }
        else
        {
            c.z = quickcos(quickinvcos(c.z) + phi1);
            c.x = quickcos(quickinvcos(c.x) + phi2);
        }
    }
    
    // the third directional cosine is fixed when normalizing
    r2.origin = r.origin;
    r2.direction = normalize(c);
    r2.diffusionAngle = 0;
    r2.source = r.source;
    return r2;
}


// shadow check.  Check if the ray from the intersection point to the light, intersection any objects inbetween.  If it does, we have a
// shadow. One improvement that can be made is to check the transparencies of all the intersecting objects and make the shadow have a
// different intensity depending on the transparency of the objects
int shadow_check(ray_s r, bubble_s *bubble, double currentTime)
{
    intersect_s i;
    ray_s r2;
    
    // to have a shadow, the light ray must intersect an object
    // first sample the ray around the shadow ray
    r2 = sample_ray(r, r.diffusionAngle);
    i = intersect(r2, bubble, currentTime);
    
    if ( i.bubble == NULL ) return 0;
    else return (r.source != i.bubble);
}


// calculate the reflection vector.
// D. Rogers: "Procedural elements for computer graphics". 5-12. p. 367
vector_s reflect(vector_s n, vector_s v1)
{
    vector_s v2;
    
    v2 = vsub(v1, svproduct(2 * vdot(n, v1), n));
    return v2;
}


// calculate the refracted vector.
// D. Rogers: "Procedural elements for computer graphics". 5-12. p. 367
vector_s refract(vector_s n, vector_s v1, double index)
{
    double p, t;
    vector_s v2;
    
    v2.x = 0.; v2.y = 0.; v2.z = 0.;
    
    p = vdot(n, v1);
    if (p < 0) {
        t = 1 - ( 1 - p*p ) / ( index*index );
        if( t <= 0 ) return v2;
        t = -p/index - sqrt(t);
    }
    else {
        index = 1 / index;
        t = 1 - ( 1 - p*p ) / ( index*index );
        if ( t <= 0 ) return v2;
        t = -p/index + sqrt(t);
    }
    v2 = vadd( svproduct(1/index, v1), svproduct(t, n) );
    return v2;
}


color_s shade(intersect_s i, ray_s r, bubble_s *bubble, double currentTime, int n);


// trace a single strait ray, n = recursion level
color_s trace_a_ray(ray_s r, bubble_s *bubble, double currentTime, int n)
{
    intersect_s inter;
    color_s col;
    double m;
    
    // check for intersection with the object
    inter = intersect(r, bubble, currentTime);
    // if no intersection, return some background color
    if (inter.bubble == NULL)
    {
        // no intersection, returning a background color
        color_s c;
        c.b = .9;
        if (r.direction.x > 0.0)
        {
            c.r = c.g = 0.4;
        }
        else
        {
            c.r = c.g = 0.4 - 0.9 * r.direction.x;
        }
        return c;
    }
    
    // calculate the shading function there
    col = shade(inter, r, bubble, currentTime, n);
    
    // if the color > 1, that means that the components are too big. Normalize.
    m = col.r;
    if (col.g > m)
        m = col.g;
    if (col.b > m)
        m = col.b;
    if (m > 1)
    {
        // overflow condition
        // normalize it!
        col.r /= m;
        col.g /= m;
        col.b /= m;
    }
    return col;
}


// trace a ray within the specified solid angle
color_s trace(ray_s r, bubble_s *bubble, double currentTime, int n)
{
    ray_s r2;
    
    r2 = sample_ray(r, r.diffusionAngle);
    return (trace_a_ray(r2, bubble, currentTime, n));
}


#define N_1 1.00029 // outside - air index of refraction
#define N_2 1.35    // soap bubble film index of refraction
#define N_3 1.00029 // inside - air index of refraction
#define R_12 (N_1 - N_2) / (N_1 + N_2)
#define R_23 (N_2 - N_3) / (N_2 + N_3)
#define EPSILON (sqr(R_12 - R_23)) / (4 * R_12 * R_23)

color_s fresnel_formula(intersect_s i, ray_s r, vector_s pointOfIntersection)
{
    color_s color;
    vector_s incident;
    double n_dot_i;
    double thickness, angle1, angle2;
    double delta, R_approx;
    double tri_X, tri_Y, tri_Z;
    double sum_x = 0, sum_y = 0, sum_z = 0, sum_Cy = 0;
    double sum_x_R = 0, sum_y_R = 0, sum_z_R = 0;
    
    thickness =  i.bubble->thicknessCoefficient * ((pointOfIntersection.y - i.bubble->yMin) / (i.bubble->yMax - i.bubble->yMin));
    incident = vsub(r.origin, pointOfIntersection);
    n_dot_i = vdot(i.normalAtPoint, incident);
    angle1 = acos(n_dot_i/(sqrt(sqr(incident.x)+sqr(incident.y)+
                                sqr(incident.z))*sqrt(sqr(i.normalAtPoint.x)+sqr(i.normalAtPoint.y)+sqr(i.normalAtPoint.z))));
    angle2 = asin((N_1 * sin(angle1)) / N_2);
    for (int j = 0, w = 380; j < 40; j++, w = w + 10) {
        delta = ((4 * M_PI) / w) * N_2 * thickness * cos(angle2);
        R_approx = (EPSILON + sqr(cos(delta / 2)));
        sum_x = sum_x + x_bar[j];
        sum_y = sum_y + y_bar[j];
        sum_z = sum_z + z_bar[j];
        sum_x_R = sum_x_R + (x_bar[j] * C[j] * R_approx);
        sum_y_R = sum_y_R + (y_bar[j] * C[j] * R_approx);
        sum_z_R = sum_z_R + (z_bar[j] * C[j] * R_approx);
        sum_Cy = sum_Cy + (y_bar[j] * C[j]);
    }
    tri_X = 4 * R_12 * R_23 * (sum_x + sum_x_R);
    tri_Y = 4 * R_12 * R_23 * (sum_y + sum_y_R);
    tri_Z = 4 * R_12 * R_23 * (sum_z + sum_z_R);
    color.r = (0.0191 * tri_X) + (-0.0053 * tri_Y) + (-0.0029 * tri_Z);
    color.g = (-0.0098 * tri_X) + (0.0200 * tri_Y) + (-0.0003 * tri_Z);
    color.b = (0.0006 * tri_X) + (-0.0012 * tri_Y) + (0.0090 * tri_Z);
    
    return color;
}


// the actual shading function.  Recursively calls trace()
color_s shade(intersect_s i, ray_s r, bubble_s *bubble, double currentTime, int n)
{
    color_s col;
    vector_s pt;
    ray_s shadow_ray;
    vector_s ldir;
    int shad;
    double ldot;
    double spec;
    vector_s rr;
    ray_s reflected; 
    color_s c;
    
    light_s light;
    light.orgin.x = -5.0;
    light.orgin.y = 8.0;
	light.orgin.z = 3.0;
    light.angle = 0.1 * M_PI/180;

    // if the recursion level has been exceeded, return peacefully
    if (n > MAX_RECURSION_LEVEL)
    {
        col.r = col.g = col.b = 0.;
        return col;
    }
    
    // initially get the ambient color
    col.r = i.bubble->ambient.r;
    col.g = i.bubble->ambient.g;
    col.b = i.bubble->ambient.b;
    
    // first calculate the intersection point
    pt = vadd(r.origin, svproduct(i.position, r.direction));
    
    color_s interferenceColor;
    interferenceColor = fresnel_formula(i, r, pt);
    
    col.r += interferenceColor.r * i.bubble->interferencePercentage;
    col.g += interferenceColor.g * i.bubble->interferencePercentage;
    col.b += interferenceColor.b * i.bubble->interferencePercentage;
    
    // for the light source first get the vector from it to the intersection
    ldir = normalize( vsub(light.orgin, pt) );
    
    // thencalc the dot product between the light direction and the normal
    // negative because the light direction is reverse (comes from the light)
    ldot = vdot(i.normalAtPoint, ldir);
    // and find if that shadow ray is stopped by an object
    shadow_ray.origin = pt;
    shadow_ray.direction = ldir;
    shadow_ray.source = i.bubble;
    // the following is the spreading angle of the ray
    shadow_ray.diffusionAngle = light.angle;
    shad = shadow_check(shadow_ray, bubble, currentTime);
    if ((ldot>0) && !shad)
    {
        // add some diffuse color
        col.r += ldot * i.bubble->diffuse.r;
        col.g += ldot * i.bubble->diffuse.g;
        col.b += ldot * i.bubble->diffuse.b;
        // now calc the specular color
        // first calculate the reflected light vector
        rr = reflect(i.normalAtPoint, ldir);
        // then take the dot of the reflected ray with the viewing direction
        // that is the specular component. The minus is there for obvious reasons
        spec = vdot(rr, r.direction);
        spec = (spec < 0) ? 0 : spec;
        // remember the specular width factor?  We use it here!
        spec = pow(spec, i.bubble->specularWidthFactor);
        col.r += spec * i.bubble->specular.r;
        col.g += spec * i.bubble->specular.g;
        col.b += spec * i.bubble->specular.b;
    }
    
    // setup the reflected ray
    reflected.origin = pt;
    reflected.direction = reflect(i.normalAtPoint, r.direction);
    reflected.source = i.bubble;
    reflected.diffusionAngle = i.bubble->reflectionDiffusion;
    
    // calculate the reflection
    c = trace(reflected, bubble, currentTime, n+1);
    
    col.r += c.r * i.bubble->reflection.r * i.bubble->reflectionPercentage;
    col.g += c.g * i.bubble->reflection.g * i.bubble->reflectionPercentage;
    col.b += c.b * i.bubble->reflection.b * i.bubble->reflectionPercentage;
        
    return col;
}


// find minimum y on bubble, n = recursion level
void find_min_y_pos(ray_s r, bubble_s *bubble, double currentTime, int n)
{
    intersect_s inter;
    vector_s point;
    
    // check for intersection with the object
    inter = intersect(r, bubble, currentTime);
    // if intersection, check if min
    if (inter.bubble != NULL)
    {
        point.y = r.direction.y * inter.position + r.origin.y;
        if (inter.bubble->yMin == 0)
            inter.bubble->yMin = point.y;
        if (inter.bubble->yMax == 0)
            inter.bubble->yMax = point.y;
        if (point.y < inter.bubble->yMin)
            inter.bubble->yMin = point.y;
        if (point.y > inter.bubble->yMax)
            inter.bubble->yMax = point.y;
    }
}


// time information
#define TIME_LIMIT_1 0.0
#define TIME_LIMIT_2 0.0

#define TRIES_PER_PIXEL 1


//  raytrace the whole scene
UIImage* raytrace(int xRes, int yRes, double bubbleThicknessCoefficient)
{
    int x, y;
    double xr, yr, xStep, yStep;
    double currentTime;
    
    color_s col, color;
    ray_s ray, r2;

    double fieldOfView = tan(37.0 * M_PI / 180) / sqrt(2.0);
    double pixelWidthInRadians = fieldOfView / xRes;
        
    vector_s up;
    up.x = 0.0;
    up.y = 1.0;
    up.z = 0.0;
    
    vector_s eye;
    eye.x = 0.0;
    eye.y = 0.0;
    eye.z = 0.4;
    
    vector_s eyeDirection;
    eyeDirection.x = 0.0;
    eyeDirection.y = 0.2;
    eyeDirection.z = 1.0;
    
    bubble_s bubble;
    bubble.center.x = 0.0;
    bubble.center.y = 1.8;
    bubble.center.z = 10.0;
    bubble.radius = 17.89;
    bubble.yMin = 0;
    bubble.yMax = 0;
    bubble.ambient.r = 0.0;
    bubble.ambient.g = 0.0;
    bubble.ambient.b = 0.0;
    bubble.diffuse.r = 0.0;
    bubble.diffuse.g = 0.0;
    bubble.diffuse.b = 0.0;
    bubble.specular.r = 1.0;
    bubble.specular.g = 1.0;
    bubble.specular.b = 1.0;
    bubble.interferencePercentage = 1.0;
    bubble.reflectionPercentage = 0.0;
    bubble.reflection.r = 0.0;
    bubble.reflection.g = 0.0;
    bubble.reflection.b = 0.0;
    bubble.specularWidthFactor = 60.0;
    bubble.reflectionDiffusion = 0.5;
    bubble.reflectionDiffusion *= M_PI / 180;
    bubble.timeMotionCoefficient.x = 0.0;
    bubble.timeMotionCoefficient.y = 0.0;
    bubble.timeMotionCoefficient.z = 0.0;
    bubble.thicknessCoefficient = bubbleThicknessCoefficient;
    NSLog(@"Coeff:%g", bubbleThicknessCoefficient);

    vector_s hor = normalize(vcross(up, eyeDirection)); // the x screen vector
    vector_s ver = normalize(vcross(eyeDirection, hor) ); // the y screen vector
    
    ray.origin = eye;  // eye is the beginning of the ray
    ray.source = NULL; // not coming from an object
    
    CGRect rect = CGRectMake(0, 0, xRes, yRes);
    UIGraphicsBeginImageContextWithOptions(rect.size, NO, 0);
    CGContextRef context = UIGraphicsGetCurrentContext();
    
    // determine the bubbles minY and maxY -- thickness increases the further down in the bubble
    yr = 1.;
    xStep = 2. / xRes; yStep = 2. / yRes;
    for(y = 0; y < yRes; y++) {
        xr = -1.;
        for(x = 0; x < xRes; x++) {
            // ray direction calculations
            ray.direction = vadd( svproduct(xr*fieldOfView, hor), svproduct(yr*fieldOfView, ver));
            ray.direction = normalize( vadd( ray.direction, eyeDirection));
            ray.diffusionAngle = 0;
            r2 = sample_ray(ray, pixelWidthInRadians);
            find_min_y_pos(r2, &bubble, 0.0, 0);
            xr += xStep;
        }
        yr -= yStep;
    }

    yr = 1.;
    xStep = 2. / xRes; yStep = 2. / yRes;
    for(y = 0; y < yRes;y++) {
        xr = -1.;
        for(x = 0; x < xRes; x++) {
            // ray direction calculations
            ray.direction = vadd(svproduct(xr*fieldOfView, hor), svproduct(yr*fieldOfView, ver));
            ray.direction = normalize(vadd(ray.direction, eyeDirection));
            ray.diffusionAngle = 0;
            col.r = col.g = col.b = 0;
            
            // the time blur has to be done in linear time and not randomly
            // randomization is used so we won't get the strobo effect
            for (int i = TRIES_PER_PIXEL; i--; ) {
                currentTime = (TIME_LIMIT_2 + TIME_LIMIT_1)/2 +
                              (TIME_LIMIT_2 - TIME_LIMIT_1) * (i + rand1()) / TRIES_PER_PIXEL;
                r2 = sample_ray(ray, pixelWidthInRadians);
                color = trace_a_ray(r2, &bubble, currentTime, 0);
                // sum all the intensities together
                col.r += color.r / TRIES_PER_PIXEL;
                col.g += color.g / TRIES_PER_PIXEL;
                col.b += color.b / TRIES_PER_PIXEL;
            }
            CGContextSetFillColorWithColor(context, [[UIColor colorWithRed:col.r green:col.g blue:col.b alpha:1.0] CGColor]);
            CGRect pixelRect = CGRectMake(x, y, x+1, y);
            CGContextFillRect(context, pixelRect);
            xr += xStep;
        }
        yr -= yStep;
    }
   
    UIImage *image = UIGraphicsGetImageFromCurrentImageContext();
    UIGraphicsEndImageContext();
    
    return image;
}
