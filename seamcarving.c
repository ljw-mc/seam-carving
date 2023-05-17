#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include "seamcarving.h"
#include "c_img.h"

int find_xlow(int x, int width){
    int xlow;

    if (x == 0){
        xlow = width - 1;
    }
    else{
        xlow = x - 1;
    }
    return xlow;    
}
int find_xhigh(int x, int width){
    int xhigh;
    if (x == width - 1){
        xhigh = 0;
    }
    else{
        xhigh = x + 1;
    }
    return xhigh;
}
int find_ylow(int y, int height){
    int ylow;
    if (y == 0){
        ylow = height - 1;
    }
    else{
        ylow = y - 1;
    }
    return ylow;
}
int find_yhigh(int y, int height){
    int yhigh;
    if (y == height - 1){
        yhigh = 0;
    }
    else{
        yhigh = y + 1;
    }
    return yhigh;

}

double sq(double x){
    return pow(x, 2);
}

uint8_t calc_grad(int rxlow, int rxhigh, int gxlow, int gxhigh, int bxlow, int bxhigh, int rylow, int ryhigh, int gylow, int gyhigh, int bylow, int byhigh){
    double Rx, Gx, Bx;
    double Ry, Gy, By;

    Rx = rxhigh - rxlow;
    Gx = gxhigh - gxlow;
    Bx = bxhigh - bxlow;

    Ry = ryhigh - rylow;
    Gy = gyhigh - gylow;
    By = byhigh - bylow;

    double sum = sq(Rx) + sq(Gx) + sq(Bx) + sq(Ry) + sq(Gy) + sq(By);
    double temp_grad = sqrt(sum);
    uint8_t grad = temp_grad/10;
    return grad;

    
}

void calc_energy(struct rgb_img *im, struct rgb_img **grad){
    // coordinates for each pixel is (y, x)
    int height = im->height;
    int width = im->width;
    int y, x, xlow, xhigh, ylow, yhigh;
    int rxlow, rxhigh, gxlow, gxhigh, bxlow, bxhigh;
    int rylow, ryhigh, gylow, gyhigh, bylow, byhigh;

    struct rgb_img ** temp_grad = malloc(3 * height * width * sizeof(struct rgb_img));
    create_img(temp_grad, height, width);

    *grad = *temp_grad;

    for (y = 0; y < height; y++){
        for (x = 0; x < width; x++) {
            /* finds the coordinates around pixel (y, x) */
            xlow = find_xlow(x, width);
            xhigh = find_xhigh(x, width);
            ylow = find_ylow(y, height);
            yhigh = find_yhigh(y, height);

            /* red */
            rxlow = get_pixel(im, y, xlow, 0);
            rxhigh = get_pixel(im, y, xhigh, 0);
            rylow = get_pixel(im, ylow, x, 0);
            ryhigh = get_pixel(im, yhigh, x, 0);
            
            /* green */
            gxlow = get_pixel(im, y, xlow, 1);
            gxhigh = get_pixel(im, y, xhigh, 1);
            gylow = get_pixel(im, ylow, x, 1);
            gyhigh = get_pixel(im, yhigh, x, 1);

            /* blue */
            bxlow = get_pixel(im, y, xlow, 2);
            bxhigh = get_pixel(im, y, xhigh, 2);
            bylow = get_pixel(im, ylow, x, 2);
            byhigh = get_pixel(im, yhigh, x, 2);

            /* calculation of grad */
            uint8_t g = calc_grad(rxlow, rxhigh, gxlow, gxhigh, bxlow, bxhigh, rylow, ryhigh, gylow, gyhigh, bylow, byhigh);
            set_pixel(*temp_grad, y, x, g, g, g);
        }
    }
}

double min3(int a, int b, int c){
    if (a <= b && a <= c) {
        return (double) a;
    } 
    else if (b <= c && b <= a) {
        return (double) b;
    } 
    else {
        return (double) c;
    }
}
double min2(int x, int y){
    if (x < y){
        return (double) x;
    }
    else{
        return (double) y;
    }
}

void dynamic_seam(struct rgb_img *grad, double **best_arr){
    int height = grad->height;
    int width = grad->width;
    int x, y, e, S;
    double *arr = malloc(3 * width * height * sizeof(double *));
    double *start = malloc(sizeof(double *));
    int center, topleft, topcenter, topright;
    start = arr;

    for (x = 0; x < width; x++){
        e = grad->raster[3*x];
        *(arr + x) = (double) e;
    }

    double temp_energy_double;
    for (y = 1; y < height; y++){
        for (x = 0; x < width; x++){
            if (x == 0){
                center = grad->raster[3*(y*width + x)];
                topcenter = *(arr + (y-1)*width + x);
                topright = *(arr + (y-1)*width + x + 1);
                temp_energy_double = min2(topcenter, topright) + (double) center;
                *(arr + y*width + x) = temp_energy_double;
            }
            else if (x == width - 1){
                center = grad->raster[3*(y*width + x)];
                topcenter = *(arr + (y-1)*width + x);
                topleft = *(arr + (y-1)*width + x - 1);
                temp_energy_double = min2(topcenter, topleft) + (double) center;
                *(arr + y*width + x) = temp_energy_double;
            }
            else{
                center = grad->raster[3*(y*width + x)];
                topcenter = *(arr + (y-1)*width + x);
                topleft = *(arr + (y-1)*width + x - 1);
                topright = *(arr + (y-1)*width + x + 1);

                temp_energy_double = min3(topleft, topcenter, topright) + (double) center;
                *(arr + y*width + x) = temp_energy_double;
            }
        }
    }

    *best_arr = start;
    return;
    
}


void recover_path(double *best, int height, int width, int **path){
    int *trail = malloc(height * sizeof(int *));
    double center, left, right;
    int min, i, k, loc;
    *path = trail;

    min = 0;
    loc = width*height - width;
    for (i = width*height - width; i < width*height; i++){
        if (*(best + i) < *(best + loc)){
            loc = i;
            min = loc - width*height + width;
        }
        *(trail + height) = min;
    }

    for (int y = height-1; y >= 0; y--){
        if (min == 0){
            center = *(best + y*width + min);
            right = *(best + y*width + (min + 1));
            if (center > right){
                min += 1;
            }
            *(trail + y) = min;
        }
        else if (min == width -1){
            center = *(best + y*width + min);
            left = *(best + y*width + (min - 1));
            if (center > left){
                min -= 1;
            }
            *(trail + y) = min;
        }
        else{
            right = *(best + y*width + (min + 1));
            center = *(best + y*width + min);
            left = *(best + y*width + (min - 1));
            if (left <= center && left <= right) {
                min -= 1;
            } 
            else if (right <= center && right <= left) {
                min += 1;
            } 
            *(trail + y) = min;
        }
    }
}

void remove_seam(struct rgb_img *src, struct rgb_img **dest, int *path){
    int x, y, reached;
    uint8_t r, g, b;
    
    int height = src->height;
    int width = src->width;

    struct rgb_img ** temp_dest = malloc(3 * height * (width-1) * sizeof(struct rgb_img *));
    create_img(temp_dest, height, width-1);
    (*temp_dest)->height = height;
    (*temp_dest)->width = width - 1;
    *dest = *temp_dest;

    // (*dest)->height = height;
    // (*dest)->width = width - 1;

    for (y = 0; y < height; y++){
        reached = 0;
        for (x = 0; x < width; x++){
            if (x == *(path + y)){
                reached = 1;
            }
            if (reached == 0 && x != *(path + y)){
                r = get_pixel(src, y, x, 0);
                g = get_pixel(src, y, x, 1);
                b = get_pixel(src, y, x, 2);

                (*temp_dest)->raster[3 * (y*((*temp_dest)->width) + x) + 0] = r;
                (*temp_dest)->raster[3 * (y*((*temp_dest)->width) + x) + 1] = g;
                (*temp_dest)->raster[3 * (y*((*temp_dest)->width) + x) + 2] = b;
            }
            else if (reached == 1 && x != *(path + y)){
                r = get_pixel(src, y, x, 0);
                g = get_pixel(src, y, x, 1);
                b = get_pixel(src, y, x, 2);

                set_pixel(*temp_dest, y, x-1, r, g, b);
            }
        }
    }
}