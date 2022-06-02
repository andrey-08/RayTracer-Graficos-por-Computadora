#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <SDL2/SDL.h>

#define WIDTH 1000
#define HEIGHT 740

#define xe 499.0
#define ye 370.0
#define ze -1000.0  

#define infinito 999999999
#define C1 0.5
#define C2 0
#define C3 0

#define x_min 0.0
#define y_min 0.0
#define EPSILON 0.0005
#define CANT_OBJ 4


typedef struct{   //Estructura utilizadas para SDL
	SDL_Renderer *render;
	SDL_Window *window;
}App;

typedef struct{   // Estructura para el Frame buffer
	double r;
	double g;
	double b;
}COLOR;

typedef struct{  //Estructura que me define un vector
	double x;
	double y;
	double z;
}PUNTO;

typedef struct{
	PUNTO centro;
	COLOR rgb;
	double R;
	double Kd; // coeficiente de reflexion, que tanto refleja la luz, proviene del archivo
	double Ka; // Coeficiente de iluminacion ambiente, proviene del archivo
	double Ks; // Propiedad de brillo especular del material del objeto
	double Kn; // Para ecualizar la reflexion especular en el material
	double O1; // que tanto color del propio objeto se observa
	double O2; // que tanto se muestra del color del objeto que se refleja
	double O3; // que tan trasparente es un objeto
	int text;
}OBJ_ESFERA;


typedef struct{
	PUNTO puntos;
	PUNTO centro;
	COLOR rgb;
	double R;
	double Kd; // coeficiente de reflexion, que tanto refleja la luz, proviene del archivo
	double Ka; // Coeficiente de iluminacion ambiente, proviene del archivo
	double Ks; // Propiedad de brillo especular del material del objeto
	double Kn; // Para ecualizar la reflexion especular en el material
	double O1; // que tanto color del propio objeto se observa
	double O2; // que tanto se muestra del color del objeto que se refleja
	double O3; // que tan trasparente es un objeto
	double D;
	int text;
}OBJ_DISCO;

typedef struct{
	PUNTO Q;
	PUNTO ancla, D2;
	COLOR rgb;
	double R;
	double Kd; // coeficiente de reflexion, que tanto refleja la luz, proviene del archivo
	double Ka; // Coeficiente de iluminacion ambiente, proviene del archivo
	double Ks; // Propiedad de brillo especular del material del objeto
	double Kn; // Para ecualizar la reflexion especular en el material
	double O1; // que tanto color del propio objeto se observa
	double O2; // que tanto se muestra del color del objeto que se refleja
	double O3; // que tan trasparente es un objeto
	int text;
}OBJ_CILINDRO;

typedef struct{
	PUNTO Q;
	PUNTO ancla, D2,D1;
	COLOR rgb;
	double R;
	double Kd; // coeficiente de reflexion, que tanto refleja la luz, proviene del archivo
	double Ka; // Coeficiente de iluminacion ambiente, proviene del archivo
	double Ks; // Propiedad de brillo especular del material del objeto
	double Kn; // Para ecualizar la reflexion especular en el material
	double O1; // que tanto color del propio objeto se observa
	double O2; // que tanto se muestra del color del objeto que se refleja
	double O3; // que tan trasparente es un objeto
	double K1, K2; // constantes para la variacion del radio segun la distancia
	int text;
}OBJ_CONOS;

typedef struct{
	PUNTO intersec;
	COLOR rgb_obj;
	PUNTO NORMAL;
	double Ka;
	double Kd;
	double Ks; 
	double Kn; 
	double t;
	double O1;
	double O2;
	double O3;
	int text;
}INTERSECTION;

typedef struct{
	PUNTO ubi;
	COLOR color;
	double Ip;
}FUENTE;

typedef struct{
	double alfa;
	double beta;
	double Y;
}CUADRATICA;

App app;
COLOR **Frame_Buffer;  	// Frame buffer
OBJ_ESFERA *esfera;    	// Listas para guardar los objetos esferas
OBJ_DISCO *discos;    	// Listas para guardar los objetos discos y planos
OBJ_CILINDRO *cilindro; // lista para guardar los objetos cillindros 
OBJ_CONOS *conos;		// Lista para guardar los objetos conos de la escena
FUENTE *luz;		   	// Lista de fuentes de luz presentes en escena

COLOR BACKGROUND;
int E_quantity, D_quantity, C_quantity, Con_quantity;
int light_quantity;
double Ia; 

void init_SDL();
void clean_window(void);
void plot_pixel(int x, int y);
void crear_buffer();
void input();
void ray_tracing();
COLOR get_color_pixel(PUNTO ancla, PUNTO direccion, int level);
INTERSECTION* get_first_intersection(PUNTO ancla, PUNTO direccion);
double angulo(double posX1, double posY1,double posZ1,double posX2, double posY2,double posZ2);
double calcular_Fatt(double dL);
void leerArchivo(int archivo);
double prod_Punto(double posX1, double posY1,double posZ1,double posX2, double posY2,double posZ2);
CUADRATICA t_conos(PUNTO ancla, PUNTO direccion, OBJ_CONOS cono);
CUADRATICA t_cilindros(PUNTO ancla, PUNTO direccion, OBJ_CILINDRO cilindro);
float frac(float pattern);
int checkerboard(int x, int y, int z);
void guardar_imagen();

int main(int argc, char *argv[]){
	int archivo = atoi(argv[1]);
	init_SDL();
	clean_window();
	crear_buffer();
	leerArchivo(archivo);
	ray_tracing();
	
	guardar_imagen();

	while (1) {    // Se va a utilizar solo para que no se cierre ventana  
		input();
    }

    SDL_DestroyWindow(app.window);
    SDL_DestroyRenderer(app.render);
    SDL_Quit();

    free(esfera);
    free(luz);
    free(cilindro);
    free(conos);
    free(discos);
	return 0;
}

void guardar_imagen(){
	FILE *img = fopen("escena1.ppm", "wb");
	(void) fprintf(img,"P6\n%d %d\n255\n",WIDTH,HEIGHT);
	for(int j = 0; j<HEIGHT; j++){
		for(int i = 0; i<WIDTH; i++){
			static unsigned char color[3];
			color[0] = Frame_Buffer[i][j].r;
			color[1] = Frame_Buffer[i][j].g;
			color[2] = Frame_Buffer[i][j].b;
			(void) fwrite(color,1,3,img);

		}
	}
	(void) fclose(img);
}

void leerArchivo(int archivo){
	FILE *escena;
	int cant_L, cant_E, cant_Con, cant_C, cant_D;
	
	if(archivo == 1){
		escena = fopen("escena1.txt", "r");
	}
	else if(archivo == 2){
		escena = fopen("escena2.txt", "r");
	}

	
	rewind(escena);

	/* Primero se extrae de la primera fila del archivo argumentos de la escena
	   como el Ia, cantidad de luces y objetos*/
	fscanf(escena, "%lf %d %d %d %d %d", &Ia, &cant_L, &cant_E, &cant_Con, &cant_C, &cant_D);

	luz = (FUENTE *)malloc(cant_L *sizeof(FUENTE));

	esfera = (OBJ_ESFERA *)malloc(cant_E * sizeof(OBJ_ESFERA));
	discos = (OBJ_DISCO *)malloc(cant_D * sizeof(OBJ_DISCO));
	cilindro = (OBJ_CILINDRO *)malloc(cant_C * sizeof(OBJ_CILINDRO));
	conos = (OBJ_CONOS *)malloc(cant_Con * sizeof(OBJ_CONOS));

	/*En las siguientes filas se obtienen primero los parametros de las fuentes de\
	  luces, como la posicion, intensidad Ip y el color*/
	for(int i = 0; i<cant_L; i++){

		double Ip,x,y,z,r,g,b;
		fscanf(escena,"%lf %lf %lf %lf %lf %lf %lf", &Ip, &x,&y,&z,&r,&g,&b);
		luz[i].Ip = Ip;
		luz[i].ubi.x = x;
		luz[i].ubi.y = y;
		luz[i].ubi.z = z;
		luz[i].color.r = r;
		luz[i].color.g = g;
		luz[i].color.b = b;
		
	}

	// Para guardar los objetos esferas //
	for(int i = 0; i<cant_E; i++){
		double ka,kd,ks,kn, O1, O2, O3,radio,xc,yc,zc,r,g,b;
		int text;
		fscanf(escena, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %i", &ka,&kd,&ks,&kn,&O1,&O2,&O3,&radio,&xc,&yc,&zc,&r,&g,&b,&text);
		
		esfera[i].centro.x=xc;
		esfera[i].centro.y=yc;
		esfera[i].centro.z=zc;
		esfera[i].rgb.r = r;
		esfera[i].rgb.g = g;
		esfera[i].rgb.b = b;
		esfera[i].R = radio;
		esfera[i].Ka = ka;
		esfera[i].Kd = kd; 
		esfera[i].Ks = ks;
		esfera[i].Kn = kn;
		esfera[i].O1 = O1;
		esfera[i].O2 = O2;
		esfera[i].O3 = O3;
		esfera[i].text = text;
	}

	// Para guardar los objetos discos y planos de la escena //
	for (int i = 0; i<cant_D; i++){
		double ka,kd,ks,kn, O1, O2, O3,radio,xc,yc,zc,r,g,b, A, B, C, D;
		int text;
		fscanf(escena, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %i", &ka,&kd,&ks,&kn,&O1,&O2,&O3,&radio,&xc,&yc,&zc,&r,&g,&b, &A,&B,&C, &D,&text);
		
		discos[i].centro.x=xc;
		discos[i].centro.y=yc;
		discos[i].centro.z=zc;
		discos[i].rgb.r = r;
		discos[i].rgb.g = g;
		discos[i].rgb.b = b;
		discos[i].R = radio;
		discos[i].Ka = ka;
		discos[i].Kd = kd; 
		discos[i].Ks = ks;
		discos[i].Kn = kn;
		discos[i].O1 = O1;
		discos[i].O2 = O2;
		discos[i].O3 = O3;
		discos[i].puntos.x= A;
		discos[i].puntos.y = B;
		discos[i].puntos.z = C;
		discos[i].D = D;
		discos[i].text = text;
	}

	// Para guardar los objetos cilindros de la escena //
	for(int i = 0; i<cant_C; i++){
		double ka,kd,ks,kn, O1, O2, O3,radio,xc,yc,zc,r,g,b, xa,ya,za, x2,y2,z2;
		int text;
		fscanf(escena, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %i", &ka,&kd,&ks,&kn,&O1,&O2,&O3,&radio,&xc,&yc,&zc,&r,&g,&b, &xa,&ya,&za, &x2,&y2,&z2,&text);
		
		cilindro[i].Q.x=xc;
		cilindro[i].Q.y=yc;
		cilindro[i].Q.z=zc;
		cilindro[i].ancla.x=xa;
		cilindro[i].ancla.y=ya;
		cilindro[i].ancla.z=za;
		cilindro[i].D2.x=x2;
		cilindro[i].D2.y=y2;
		cilindro[i].D2.z=z2;
		cilindro[i].rgb.r = r;
		cilindro[i].rgb.g = g;
		cilindro[i].rgb.b = b;
		cilindro[i].R = radio;
		cilindro[i].Ka = ka;
		cilindro[i].Kd = kd; 
		cilindro[i].Ks = ks;
		cilindro[i].Kn = kn;
		cilindro[i].O1 = O1;
		cilindro[i].O2 = O2;
		cilindro[i].O3 = O3;
		cilindro[i].text = text;
	}

	//Para guardar los objetos conos de la escena //
	for(int i = 0; i<cant_Con; i++){
		double ka,kd,ks,kn, O1, O2, O3,radio,xc,yc,zc,r,g,b, xa,ya,za, x2,y2,z2, x1,y1,z1,k1,k2;
		int text;
		fscanf(escena, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %i", &ka,&kd,&ks,&kn,&O1,&O2,&O3,&radio,&xc,&yc,&zc,&r,&g,&b, &xa,&ya,&za,&x2,&y2,&z2,&x1,&y1,&z1,&k1,&k2,&text);
		
		conos[i].Q.x=xc;
		conos[i].Q.y=yc;
		conos[i].Q.z=zc;
		conos[i].ancla.x=xa;
		conos[i].ancla.y=ya;
		conos[i].ancla.z=za;
		conos[i].D2.x=x2;
		conos[i].D2.y=y2;
		conos[i].D2.z=z2;
		conos[i].D1.x=x1;
		conos[i].D1.y=y1;
		conos[i].D1.z=z1;
		conos[i].rgb.r = r;
		conos[i].rgb.g = g;
		conos[i].rgb.b = b;
		conos[i].R = radio;
		conos[i].Ka = ka;
		conos[i].Kd = kd; 
		conos[i].Ks = ks;
		conos[i].Kn = kn;
		conos[i].O1 = O1;
		conos[i].O2 = O2;
		conos[i].O3 = O3;
		conos[i].K1 = k1;
		conos[i].K2 = k2;
		conos[i].text = text;
	}

	E_quantity = cant_E;
	D_quantity = cant_D;
	C_quantity = cant_C;
	Con_quantity = cant_Con;
	light_quantity = cant_L;

	fclose(escena);
}

void ray_tracing(){ 
	COLOR color;
	PUNTO direc;
	PUNTO ancla;
	ancla.x = xe;
	ancla.y = ye;
	ancla.z = ze;

	//while (1){
	    for(int i = 0; i < WIDTH; i++){ //for para el frame buffer hacia abajo
	        for(int j = 0; j < HEIGHT; j++){//for para el frame buffer hacia derecha 
	            double dx = (double) i; 
	            double dy = (double) j;
	            double xw = (double)((dx + 0.5)*(WIDTH - x_min) / WIDTH + x_min);
	            double yw = (double)((dy + 0.5)*(HEIGHT - y_min) / HEIGHT + y_min);
	            double zw = 0.0; 
	 
	            double mag = sqrt(pow((xw-ancla.x),2)+pow((yw-ancla.y),2)+pow((zw-ancla.z),2));
	            direc.x = (xw - ancla.x) / mag;
	            direc.y = (yw - ancla.y) / mag;
	            direc.z = (zw - ancla.z) / mag;
	            color = get_color_pixel(ancla,direc, 6);

	            Frame_Buffer[i][j].r = color.r; 
	            Frame_Buffer[i][j].g = color.g; 
	            Frame_Buffer[i][j].b = color.b;

	            SDL_SetRenderDrawColor(app.render, Frame_Buffer[i][j].r, Frame_Buffer[i][j].g,Frame_Buffer[i][j].b, 0);
	            plot_pixel(i,j);
	            
	            // SDL_RenderPresent(app.render);
	            // input();
	        }
	    }
	    SDL_RenderPresent(app.render);
	   // exit(0);
	//}
}

int checkerboard(int x, int y, int z){
	float pattern = (x+y);
	pattern = frac(pattern/2);
	pattern *= 2;
	return (int)pattern;
}

float frac(float pattern){
	double frac, entera;
	frac = modf(pattern, &entera);
	return frac;
}

COLOR get_color_pixel(PUNTO ancla, PUNTO direccion, int level){
	COLOR color;
	INTERSECTION *interseccion;
	interseccion = get_first_intersection(ancla,direccion);
	INTERSECTION *obstaculo; 
	PUNTO N;
	PUNTO L;
	PUNTO V;
	PUNTO R;
	PUNTO reflejo;

	COLOR color_reflejado, color_transp;

	if(!interseccion){
		color.r = 50;
		color.g = 50;
		color.b = 50;
	}
	else{

		double I = 0;
		double E = 0;

		double inter_x = interseccion->intersec.x;
		double inter_y = interseccion->intersec.y;
		double inter_z = interseccion->intersec.z;

		N.x = interseccion -> NORMAL.x;
		N.y = interseccion -> NORMAL.y;
		N.z = interseccion -> NORMAL.z;
		
		//Calculamos el vector V (normalizado) //
		V.x = -1*(direccion.x);
		V.y = -1*(direccion.y);
		V.z = -1*(direccion.z);	

		for(int k = 0; k<light_quantity; k++){
			// Averiguando el vector L //
			double norma_L = sqrt(pow((luz[k].ubi.x-inter_x),2)+pow((luz[k].ubi.y-inter_y),2)+pow((luz[k].ubi.z-inter_z),2));
			L.x=(luz[k].ubi.x - inter_x)/norma_L;
			L.y=(luz[k].ubi.y - inter_y)/norma_L;
			L.z=(luz[k].ubi.z - inter_z)/norma_L;
			double productoPunto = angulo(L.x,L.y,L.z,N.x,N.y,N.z);

			// Averiguando el vector R //
			R.x = 2*(N.x)*(productoPunto) - (L.x);
			R.y = 2*(N.y)*(productoPunto) - (L.y);
			R.z = 2*(N.z)*(productoPunto) - (L.z);	
			double productoEspecular =  angulo(R.x,R.y,R.z,V.x,V.y,V.z);
 
			if(productoPunto > 0 ){ //Toma en cuenta la fuente de luz
				double Fatt = calcular_Fatt(norma_L);
				obstaculo = get_first_intersection(interseccion->intersec, L);

				if(!obstaculo || obstaculo->t > norma_L ){
					I = I + (productoPunto*(interseccion->Kd)*Fatt*(luz[k].Ip));

					if(productoEspecular > 0){
						E = E + (pow(productoEspecular,interseccion->Kn)*(interseccion->Ks)*Fatt*(luz[k].Ip));
					}
				}
				free(obstaculo);
			}
		}
		
		I = I + (Ia*(interseccion->Ka)); // Intensidad final con respecto a la fuente de luz
		if(I>1){
			I=1;
		}

		if(E>1){  // Forzamos a E para que sea 1 en caso de que sobrepase 1
			E=1;
		}

		if(interseccion -> text == 1){
			int num = checkerboard((int)inter_x,(int)inter_y,(int)inter_z);
			double rojo, verde, azul;
			if(num % 2 == 0){ // si es par pinta blanco
				rojo = 255;
				verde = 0;
				azul = 0;
			}
			else{
				rojo = 0;
				verde = 0;
				azul = 0;
			} 
			color.r = I*rojo;
			color.g = I*verde;
			color.b = I*azul;
		}

		else{
			color.r = I*(interseccion ->rgb_obj.r);
			color.g = I*(interseccion ->rgb_obj.g);
			color.b = I*(interseccion ->rgb_obj.b);
		}

		// AGREGAMOS REFLEXION ESPECULAR  //
		color.r = color.r + E*(255 - color.r);
		color.g = color.g + E*(255 - color.g);
		color.b = color.b + E*(255 - color.b);

		/* AQUI SE CALCULA EL COLOR REFLEJADO */
		if(interseccion -> O2 != 0){
			double prod_ND = angulo(N.x,N.y,N.z,direccion.x,direccion.y,direccion.z);
			reflejo.x = -2*(N.x)*(prod_ND) - (direccion.x);
			reflejo.y = -2*(N.y)*(prod_ND) - (direccion.y);
			reflejo.z = -2*(N.z)*(prod_ND) - (direccion.z);

			if(level > 0){
				color_reflejado = get_color_pixel(interseccion->intersec,reflejo, level - 1);
			}
		}

		// AQUI IRIA TRANSPARENCIA //

		if(interseccion -> O3 != 0){
			if(level>0){
				color_transp = get_color_pixel(interseccion->intersec, direccion, level - 1);
			}
		}


		color.r = (interseccion -> O1)*color.r + (interseccion->O2)*color_reflejado.r + (interseccion->O3)*color_transp.r;
		color.g = (interseccion -> O1)*color.g + (interseccion->O2)*color_reflejado.g + (interseccion->O3)*color_transp.g; 
		color.b = (interseccion -> O1)*color.b + (interseccion->O2)*color_reflejado.b + (interseccion->O3)*color_transp.b;

		free(interseccion);
	}

	return (color); // Funcion para obtener el color del pixel a pintar
}

double calcular_Fatt(double dL){
	double Fatt;
	Fatt = 1/(C1 + (C2*dL) + (C3*pow(dL,2)));
	if(0 <= Fatt && Fatt<=1){
		return (Fatt);
	}
	else if(Fatt > 1){
		return (1.0);
	}
}

double angulo(double posX1, double posY1,double posZ1,double posX2, double posY2,double posZ2){
	double productoPunto = (posX1*posX2) + (posY1*posY2) + (posZ2*posZ1);
	double magnitud1 = sqrt(pow(posX1,2) + pow(posY1,2) + pow(posZ1,2));
    double magnitud2 = sqrt(pow(posX2,2) + pow(posY2,2) + pow(posZ2,2)); 
    return (productoPunto/(magnitud1*magnitud2));
}

double prod_Punto(double posX1, double posY1,double posZ1,double posX2, double posY2,double posZ2){
	return (posX1*posX2 + posY1*posY2 + posZ1*posZ2);
}

CUADRATICA t_conos(PUNTO ancla, PUNTO direccion, OBJ_CONOS cono){
	double xq, yq, zq, x0 , y0, z0;
	CUADRATICA cuadratica;
	xq = cono.Q.x;
	yq = cono.Q.y; 
	zq = cono.Q.z;

	x0 = cono.ancla.x;
	y0 = cono.ancla.y;
	z0 = cono.ancla.z;

	double DQ = direccion.x * xq + direccion.y*yq + direccion.z*zq;
	double AQ = ancla.x*xq + ancla.y*yq + ancla.z*zq;
	double CQ = x0*xq + y0*yq + z0*zq;

	double A = xq*DQ - direccion.x;
	double B = x0 - ancla.x + xq*(AQ - CQ);
	double C = yq*DQ - direccion.y;
	double D = y0 - ancla.y + yq*(AQ - CQ);
	double E = zq*DQ - direccion.z;
	double F = z0 - ancla.z + zq*(AQ - CQ);
	double G = (cono.K2/cono.K1)*DQ;
	double H = (cono.K2/cono.K1)*(AQ - CQ);

	cuadratica.alfa = A*A + C*C + E*E - G*G;
	cuadratica.beta = 2*(A*B + C*D + E*F - G*H);
	cuadratica.Y = B*B + D*D + F*F - H*H;

	return cuadratica;
}

CUADRATICA t_cilindros(PUNTO ancla, PUNTO direccion, OBJ_CILINDRO cilindro){
	double xq, yq, zq, x0 , y0, z0;
	CUADRATICA cuadratica;
	xq = cilindro.Q.x;
	yq = cilindro.Q.y; 
	zq = cilindro.Q.z;

	x0 = cilindro.ancla.x;
	y0 = cilindro.ancla.y;
	z0 = cilindro.ancla.z;

	double DQ = direccion.x * xq + direccion.y*yq + direccion.z*zq;
	double AQ = ancla.x*xq + ancla.y*yq + ancla.z*zq;
	double CQ = x0*xq + y0*yq + z0*zq;

	double A = xq*DQ - direccion.x;
	double B = x0 - ancla.x + xq*(AQ - CQ);
	double C = yq*DQ - direccion.y;
	double D = y0 - ancla.y + yq*(AQ - CQ);
	double E = zq*DQ - direccion.z;
	double F = z0 - ancla.z + zq*(AQ - CQ);

	cuadratica.alfa = A*A + C*C + E*E;
	cuadratica.beta = 2*(A*B + C*D + E*F);
	cuadratica.Y = B*B + D*D + F*F - cilindro.R * cilindro.R;

	return cuadratica;
}

INTERSECTION* get_first_intersection(PUNTO ancla, PUNTO direccion){

	INTERSECTION *interseccion = NULL; // Valor a retornar

	PUNTO int_point;
	COLOR rgb;
	PUNTO centro_int;
	PUNTO normal;

	long double t_min = infinito;
	double xc,yc,zc; //Coordaenadas del centro de la esfera y vector Q del cilindro
	double xa,ya,za; 
	double radio, radioF;    //Radio de la esfera
	double beta, Y, alfa;  //Para calcular discriminante e intersecciones
	double discrim;  //discriminante
	int detected = 0; // Para informar cuando detectamos una interseccion
	double k_a,k_d,k_s,k_n;
	double o1,o2, o3;
	int txt;
	int cant;

	for(int i = 0; i<CANT_OBJ; i++){ // de tipo a tipo de objetos

		if(i == 0){
			cant = E_quantity;
		}
		if( i == 1){
			cant = D_quantity;
		}
		if( i == 2){
			cant = C_quantity;
		}

		if( i == 3){
			cant = Con_quantity;
		}

		for(int j = 0; j<cant; j++){


			double t1 = infinito;
			double t2 = infinito;
			double distancia = infinito;
		    
			if( i == 0 ){ // para una esfera
				xc = esfera[j].centro.x;
				yc = esfera[j].centro.y;
				zc = esfera[j].centro.z;
				radio = esfera[j].R;	
				beta = 2*(direccion.x * (ancla.x - xc) + direccion.y * (ancla.y - yc) + direccion.z * (ancla.z - zc));
				Y = pow((ancla.x - xc), 2) + pow((ancla.y - yc),2) + pow((ancla.z - zc), 2) - pow(radio, 2);
				alfa = 1;
				discrim = pow(beta, 2) - (4.0 * Y*alfa);
			}

			else if( i == 1){  // Para discos y planos
				double num = (discos[j].puntos.x*ancla.x + discos[j].puntos.y*ancla.y + discos[j].puntos.z*ancla.z + discos[j].D);
				double denom = (discos[j].puntos.x*direccion.x + discos[j].puntos.y*direccion.y + discos[j].puntos.z*direccion.z);
				t1 = -num/denom;
				if(t1 > 0){
					distancia = t1;
				}
			}

			else if(i == 2){ // Para un cilindro

				radio = cilindro[j].R;
				CUADRATICA aux1;
				aux1 = t_cilindros(ancla , direccion, cilindro[j]);
				alfa = aux1.alfa;
				beta = aux1.beta;
				Y = aux1.Y;
				discrim = pow(beta,2) - (4*alfa*Y);

			}

			else if (i == 3){  // Para los conos
				CUADRATICA aux;
				aux = t_conos(ancla , direccion, conos[j]);
				alfa = aux.alfa;
				beta = aux.beta;
				Y = aux.Y;

				discrim = pow(beta,2) - (4*alfa*Y);		
			}
			

			if(i != 1 ){	// Para esfera y cilindro
				if(discrim >= 0){
					if(discrim > EPSILON){
						t1 = ((-beta + sqrt(pow(beta,2) - 4 * Y* alfa))/(2*alfa));
						t2 = ((-beta - sqrt(pow(beta,2) - 4 * Y * alfa))/(2*alfa));
						// printf("t1: %f  t2: %f \n", t1, t2);

						if(t1>t2){
							distancia = t2;
						}
						else if(t2>t1){
							distancia = t1;
						}
					}
					else if (discrim >= 0 && discrim <= EPSILON ){
						t1 = ((-beta + sqrt(pow(beta,2) - 4 * Y))/2);
						distancia = t1;
					}
				}
			}

			// Entre todos los objetos decide cual es el que tiene la interseccion de menor distancia //
			if(distancia < t_min  && distancia>0 && distancia>EPSILON){
				double xi, yi,zi;

				// Guardamos parametros de la esfera //
				if( i== 0){
					// Calculamos las intersecciones con el objeto //
					int_point.x = ancla.x + distancia * direccion.x;
					int_point.y = ancla.y + distancia * direccion.y;
					int_point.z = ancla.z + distancia * direccion.z;

					centro_int.x=xc;
					centro_int.y=yc;
					centro_int.z=zc;
					radioF = radio;

					rgb.r = esfera[j].rgb.r;
					rgb.g = esfera[j].rgb.g;
					rgb.b = esfera[j].rgb.b;
					k_a = esfera[j].Ka;
					k_d = esfera[j].Kd;
					k_n = esfera[j].Kn;
					k_s = esfera[j].Ks;
					o1 = esfera[j].O1;
					o2 = esfera[j].O2;
					o3 = esfera[j].O3;
					txt = esfera[j].text;


					double norma_N = radioF;
					normal.x= (int_point.x - centro_int.x)/norma_N;
					normal.y= (int_point.y - centro_int.y)/norma_N;
					normal.z= (int_point.z - centro_int.z)/norma_N;
					t_min = distancia;
				}

				// Guardamos parametros del disco //
				else if( i == 1){

					// Calculamos las intersecciones con el objeto //
					xi = ancla.x + distancia * direccion.x;
					yi = ancla.y + distancia * direccion.y;
					zi = ancla.z + distancia * direccion.z;
					double delimited = sqrt(pow(xi - discos[j].centro.x,2)+pow(yi - discos[j].centro.y,2)+pow(zi - discos[j].centro.z,2));

					if(delimited <= discos[j].R){
						// Calculamos las intersecciones con el objeto //
						int_point.x = xi;
						int_point.y = yi;
						int_point.z = zi;

						rgb.r = discos[j].rgb.r;
						rgb.g = discos[j].rgb.g;
						rgb.b = discos[j].rgb.b;
						k_a = discos[j].Ka;
						k_d = discos[j].Kd;
						k_n = discos[j].Kn;
						k_s = discos[j].Ks;
						o1 = discos[j].O1;
						o2 = discos[j].O2;
						o3 = discos[j].O3;
						txt = discos[j].text;
						normal = discos[j].puntos;
						t_min = distancia;
					}
				}

				// Guardamos parametros del cilindro //
				else if(i == 2){
					radioF = radio;
					
					// Calculamos las intersecciones con el objeto //
					xi = ancla.x + distancia * direccion.x;
					yi = ancla.y + distancia * direccion.y;
					zi = ancla.z + distancia * direccion.z;
					double d = (xi - cilindro[j].ancla.x)*cilindro[j].Q.x + (yi - cilindro[j].ancla.y)*cilindro[j].Q.y + (zi - cilindro[j].ancla.z)*cilindro[j].Q.z;
					double delimited = sqrt(pow(cilindro[j].D2.x - cilindro[j].ancla.x,2)+pow(cilindro[j].D2.y - cilindro[j].ancla.y,2)+pow(cilindro[j].D2.z - cilindro[j].ancla.z,2));
					
					if(d <= delimited && d>=0){
						// Calculamos las intersecciones con el objeto //
						
						int_point.x = xi;
						int_point.y = yi;
						int_point.z = zi;

						rgb.r = cilindro[j].rgb.r;
						rgb.g = cilindro[j].rgb.g;
						rgb.b = cilindro[j].rgb.b;
						k_a = cilindro[j].Ka;
						k_d = cilindro[j].Kd;
						k_n = cilindro[j].Kn;
						k_s = cilindro[j].Ks;
						o1 = cilindro[j].O1;
						o2 = cilindro[j].O2;
						o3 = cilindro[j].O3;
						txt = cilindro[j].text;

						double xm,ym,zm;
						// calculando la normal //
						xm = cilindro[j].ancla.x + d*cilindro[j].Q.x; 
						ym = cilindro[j].ancla.y + d*cilindro[j].Q.y;
						zm = cilindro[j].ancla.z + d*cilindro[j].Q.z;

						normal.x= (int_point.x - xm)/radioF;
						normal.y= (int_point.y - ym)/radioF;
						normal.z= (int_point.z - zm)/radioF;
						t_min = distancia;
					}

				}

				else if( i == 3){// Guardamos parametros del cono //
					double xm,ym,zm;

					// Calculamos las intersecciones con el objeto //
					xi = ancla.x + distancia * direccion.x;
					yi = ancla.y + distancia * direccion.y;
					zi = ancla.z + distancia * direccion.z;
					
					double delimited = sqrt(pow(conos[j].ancla.x - conos[j].D1.x,2)+pow(conos[j].ancla.y - conos[j].D1.y,2)+pow(conos[j].ancla.z - conos[j].D1.z,2));
					double delimited2 = sqrt(pow(conos[j].ancla.x - conos[j].D2.x,2)+pow(conos[j].ancla.y - conos[j].D2.y,2)+pow(conos[j].ancla.z - conos[j].D2.z,2));
					double d_ratio = (xi-conos[j].ancla.x)*conos[j].Q.x + (yi-conos[j].ancla.y)*conos[j].Q.y +(zi-conos[j].ancla.z)*conos[j].Q.z;
					
					if(d_ratio >= delimited && d_ratio<= delimited2){
						int_point.x = xi;
						int_point.y = yi;
						int_point.z = zi;

						rgb.r = conos[j].rgb.r;
						rgb.g = conos[j].rgb.g;
						rgb.b = conos[j].rgb.b;
						k_a = conos[j].Ka;
						k_d = conos[j].Kd;
						k_n = conos[j].Kn;
						k_s = conos[j].Ks;
						o1 = conos[j].O1;
						o2 = conos[j].O2;
						o3 = conos[j].O3;
						txt = conos[j].text;

						// calculando la normal //
						double xm,ym,zm, mag_n;
						double mag_L = sqrt(pow(xi - conos[j].ancla.x,2)+pow(yi - conos[j].ancla.y,2)+pow(zi - conos[j].ancla.z,2));
						double d = (mag_L*mag_L) / ((xi - conos[j].ancla.x)*conos[j].Q.x +(yi - conos[j].ancla.y)*conos[j].Q.y +(zi - conos[j].ancla.z)*conos[j].Q.z);
						xm = conos[j].ancla.x + d*conos[j].Q.x; 
						ym = conos[j].ancla.y + d*conos[j].Q.y;
						zm = conos[j].ancla.z + d*conos[j].Q.z;

						mag_n = sqrt(pow(int_point.x - xm,2)+ pow(int_point.y - ym,2)+pow(int_point.z - zm,2));
						normal.x= (int_point.x - xm)/mag_n;
						normal.y= (int_point.y - ym)/mag_n;
						normal.z= (int_point.z - zm)/mag_n;
						t_min = distancia;
					}
				}	

				detected = 1;
			}	
		}	
	}

	if(detected){
		//printf("Hola\n");
		interseccion= (INTERSECTION *)malloc(sizeof(INTERSECTION));
		interseccion -> intersec = int_point; // Aqui se iguala a la interseccion calculada
		interseccion -> rgb_obj = rgb;
		interseccion -> Ka = k_a;
		interseccion -> Kd = k_d;
		interseccion -> Kn = k_n;
		interseccion -> Ks = k_s;
		interseccion -> t = t_min;
		interseccion -> O1 = o1;
		interseccion -> O2 = o2;
		interseccion -> O3 = o3;
		interseccion -> text = txt;
		interseccion -> NORMAL = normal;

	}

	return (interseccion);
}

void crear_buffer(){  // Funcion par ainicializar y reservar espacio para frame buffer
	
	if(Frame_Buffer != NULL){  // Si el buffer ya esta inicializado

		for(int i = 0;  i< WIDTH; i++){ // Se liberan las filas que es donde estan contenidas las columnas
			free(Frame_Buffer[i]);
		}
		free(Frame_Buffer);
	}

	// Reservamos espacio en memoria para el frame buffer //
	Frame_Buffer = (COLOR **)malloc(WIDTH * sizeof(COLOR *));
	for(int i= 0 ; i<WIDTH; i++){
		Frame_Buffer[i]= (COLOR *)malloc(HEIGHT * sizeof(COLOR));
	}

	// Inicializamos el frame con color blanco en todos los pixeles //
	for( int i = 0; i<WIDTH; i++){
		for(int j = 0; j<HEIGHT; j++){
			Frame_Buffer[i][j].r = 255;
			Frame_Buffer[i][j].g = 255;
			Frame_Buffer[i][j].b = 255;
		}
	}
}

void input(void){ // Para gestionar entradas o salidas de SDL
	SDL_Event event;

	while(SDL_PollEvent(&event)){		
		switch(event.type){
			case SDL_QUIT:
				exit(0);
				break;
			default:
				break;
		}

	}
}

void init_SDL(){
	SDL_Init(SDL_INIT_VIDEO);
	app.window = SDL_CreateWindow("Ray Tracer / ITCR",  SDL_WINDOWPOS_CENTERED,  SDL_WINDOWPOS_CENTERED,WIDTH, HEIGHT, 0);
	app.render = SDL_CreateRenderer(app.window, -1, SDL_RENDERER_ACCELERATED);
}

void clean_window(void){
	SDL_SetRenderDrawColor(app.render, 0, 0, 0, 0);
	SDL_RenderClear(app.render);
	SDL_RenderPresent(app.render); // Para limpiar la ventana
}

void plot_pixel(int x, int y){	
	SDL_RenderDrawPoint(app.render, x, y);
}