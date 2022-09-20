#define _CRT_SECURE_NO_WARNINGS
#define DOKLADNOSC 31 // nieparzyste! jakie k i l maksymalne maj¹ byæ
#define M_PI 3.14159265358979323846
#include <stdio.h>
#include <math.h>
#include <locale.h>
//#include <time.h>

double lepszyScanfDoDoublea()
{
	double liczba;
	int c;
	int i = 0;

	while (1)
	{
		if (!scanf("%lf", &liczba))
		{
			while ('\n' != getchar());
			printf("Sprobuj ponownie: ");
		}
		else
		{
			while ((c = getchar()) != '\n' && c != EOF)
			{
				i++;
			}
			if (i > 0)
			{
				printf("Wskazówka: jako separatora dziesiêtnego u¿ywaj przecinka.\nSprobuj ponownie: ");
				liczba = lepszyScanfDoDoublea();
			}
			return liczba;
		}


	}
}

int lepszyScanfDoInta()
{
	double liczbaDowolna = 0;
	int liczbaCalkowita = 0;
	int i = 1;
	while (i)
	{
		liczbaDowolna = lepszyScanfDoDoublea();
		if (liczbaDowolna == (floor(liczbaDowolna)))
		{
			liczbaCalkowita = (int)liczbaDowolna;
			i = 0;
		}
		else printf("\nWymagana liczba calkowita!!!\nSprobuj ponownie: ");
	}
	return liczbaCalkowita;
}

void wypiszKolor(FILE* plik, double naprezenie, double naprezenieMax, int glebia)
{
	int czerwony, zielony, niebieski;
	int identyfikator = round(naprezenie / naprezenieMax * 4 * glebia);
	if (identyfikator < glebia)
	{
		czerwony = 0;
		zielony = identyfikator;
		niebieski = glebia;
	}
	else if (identyfikator < 2 * glebia)
	{
		czerwony = 0;
		zielony = glebia;
		niebieski = 2 * glebia - identyfikator;
	}
	else if (identyfikator < 3 * glebia)
	{
		czerwony = identyfikator - 2 * glebia;
		zielony = glebia;
		niebieski = 0;
	}
	else
	{
		czerwony = glebia;
		zielony = 4 * glebia - identyfikator;
		niebieski = 0;
	}
	fprintf(plik, "%d %d %d", czerwony, zielony, niebieski);
}

void tauObliczanie(double* tau, double a, double b, double Ms, double J, double xm, double ym)
{
	int k, l;
	double tymczasowa, b2, sl[(DOKLADNOSC + 1) / 2], cl[(DOKLADNOSC + 1) / 2], l2a2[(DOKLADNOSC + 1) / 2], b2l[(DOKLADNOSC + 1) / 2], a2l3[(DOKLADNOSC + 1) / 2];
	b2 = b * b;
	tau[0] = 0; tau[1] = 0;
	if (a >= b)
	{
		for (l = 1; l <= (DOKLADNOSC + 1); l = l + 2)
		{
			sl[(l - 1)/2] = sin(l * M_PI * ym / b);
			cl[(l - 1) / 2] = cos(l * M_PI * ym / b);
			l2a2[(l - 1) / 2] = l * l * a * a;
			b2l[(l - 1) / 2] = b * b * l;
			a2l3[(l - 1) / 2] = a * a * l * l * l;
		}
		for (k = 1; k <= DOKLADNOSC; k = k + 2)
		{
			for (l = 0; l < (DOKLADNOSC + 1) / 2; l = l + 1)
			{
				tau[0] = tau[0] + pow(-1, (k + 2 * l + 1) / 2)*cos(k*M_PI*xm / a)*sl[l] / (k*k*k*b2 + l2a2[l] * k);
				tau[1] = tau[1] + pow(-1, (k + 2 * l + 1) / 2)*sin(k*M_PI*xm / a)*cl[l] / (k*k*b2l[l] + a2l3[l]);
			}
		}
		tau[0] = Ms * pow(2, 5) * a * a * b / J / pow(M_PI, 3) * tau[0];
		tau[1] = -Ms * pow(2, 5) * a * b * b / J / pow(M_PI, 3) * tau[1];
	}
	else
	{
		tymczasowa = a; a = b; b = tymczasowa;
		tymczasowa = xm; xm = ym; ym = -tymczasowa;
		for (l = 1; l <= (DOKLADNOSC + 1); l = l + 2)
		{
			sl[(l - 1) / 2] = sin(l * M_PI * ym / b);
			cl[(l - 1) / 2] = cos(l * M_PI * ym / b);
			l2a2[(l - 1) / 2] = l * l * a * a;
			b2l[(l - 1) / 2] = b * b * l;
			a2l3[(l - 1) / 2] = a * a * l * l * l;
		}
		for (k = 1; k <= DOKLADNOSC; k = k + 2)
		{
			for (l = 0; l < (DOKLADNOSC + 1) / 2; l = l + 1)
			{
				tau[0] = tau[0] + pow(-1, (k + 2 * l + 1) / 2)*cos(k*M_PI*xm / a)*sl[l] / (k*k*k*b2 + l2a2[l] * k);
				tau[1] = tau[1] + pow(-1, (k + 2 * l + 1) / 2)*sin(k*M_PI*xm / a)*cl[l] / (k*k*b2l[l] + a2l3[l]);
			}
		}
		tau[0] = Ms * pow(2, 5) * a * a * b / J / pow(M_PI, 3) * tau[0];
		tau[1] = -Ms * pow(2, 5) * a * b * b / J / pow(M_PI, 3) * tau[1];
		tymczasowa = tau[0]; tau[0] = -tau[1]; tau[1] = tymczasowa;
	}
}

double K1(double r)
{
	double suma = 0;
	double l4r2[(DOKLADNOSC + 1) / 2], r2;
	int k, l, l2[(DOKLADNOSC + 1) / 2];
	r2 = r * r;
	for (l = 1; l <= (DOKLADNOSC + 1) / 2; l = l + 1)
	{
		l2[l - 1] = l * l;
		l4r2[l - 1] = l * l * l * l * r * r;
	}
	for (k = 1; k <= DOKLADNOSC; k = k + 2)
	{
		for (l = 0; l < (DOKLADNOSC + 1) / 2; l = l + 1)
		{
			suma = suma + r2 / (l2[l] * pow(k, 4) + pow(k, 2)*l4r2[l]);
		}
	}
	return suma * pow(2, 8) / pow(M_PI, 6);
}

void wczytywaniePodstawowychDanych(double* Ms, double* Mgx, double* Mgy, double* T, double* P, int* wymiar, int* glebia)
{
	printf("Wczytywanie danych\n\nObci¹¿enie\nMs[Nm] = "); *Ms = lepszyScanfDoDoublea(); // os zwrocona do uzytkownika
	printf("Mgx[Nm]= ");	*Mgx = lepszyScanfDoDoublea(); // os zwrocona w prawo
	printf("Mgy[Nm]= ");	*Mgy = lepszyScanfDoDoublea(); // os zwrocona do gory
	printf("Ty[N]  = ");	*T = lepszyScanfDoDoublea(); // dodatnie w dol
	printf("Pz[N]  = ");	*P = lepszyScanfDoDoublea(); // dodatnie rozciaganie
	printf("\nParametry grafiki:\t");
	while (1)
	{
		printf("\nWymiar obrazka: "); 	*wymiar = lepszyScanfDoInta();
		if (*wymiar > 0) break;
		else printf("Wymiar musi byæ dodatni!\n");
	}
	while (1)
	{
		printf("G³êbia koloru:  "); 	*glebia = lepszyScanfDoInta();
		if (*glebia > 0) break;
		else printf("G³êbia musi byæ dodatnia!\n");
	}
	(*glebia)--;
	printf("\nWymiary geometryczne: \n");
}

void przekrojKolowy()
{
	FILE* plik;
	double Ms, Mgx, Mgy, d, T, P;
	double naprezenieMax = 0;
	double naprezenie;
	int x, y, wymiar, glebia;
	wczytywaniePodstawowychDanych(&Ms, &Mgx, &Mgy, &T, &P, &wymiar, &glebia);
	while (1)
	{
		printf("d[m]: ");	scanf("%lf", &d);
		if (d > 0) break;
		else printf("Œrednica musi byæ dodatnia!\n");
	}
	for (y = 0; y < wymiar; y++)
	{
		for (x = 0; x < wymiar; x++)
		{
			if ((pow(wymiar / 2 - x, 2) + pow(wymiar / 2 - y, 2)) <= (wymiar * wymiar / 4))
			{
				naprezenie = 32 / M_PI / pow(d, 4)*sqrt(pow(2 * (Mgx*(wymiar / 2 - y)*d / wymiar - Mgy * (x - wymiar / 2)*d / wymiar) + P * pow(d, 2) / 8, 2) + 3 * (pow(Ms*(x - wymiar / 2)*d / wymiar, 2) + pow(-Ms * (wymiar / 2 - y)*d / wymiar - T / 3 * (d*d / 4 - (wymiar / 2 - y)*(wymiar / 2 - y)*d*d / wymiar / wymiar), 2)));
				if (naprezenie > naprezenieMax) naprezenieMax = naprezenie;
			}
		}
	}
	plik = fopen("wizualizacja - przekrój ko³owy.ppm", "w");
	if (plik != NULL)
	{
		fprintf(plik, "P3\n%d %d\n%d\n", wymiar, wymiar, glebia);
		for (y = 0; y < wymiar; y++)
		{
			for (x = 0; x < wymiar; x++)
			{
				if ((pow(wymiar / 2 - x, 2) + pow(wymiar / 2 - y, 2)) <= (wymiar * wymiar / 4))
				{
					naprezenie = 32 / M_PI / pow(d, 4)*sqrt(pow(2 * (Mgx*(wymiar / 2 - y)*d / wymiar - Mgy * (x - wymiar / 2)*d / wymiar) + P * pow(d, 2) / 8, 2) + 3 * (pow(Ms*(x - wymiar / 2)*d / wymiar, 2) + pow(-Ms * (wymiar / 2 - y)*d / wymiar - T / 3 * (d*d / 4 - (wymiar / 2 - y)*(wymiar / 2 - y)*d*d / wymiar / wymiar), 2)));
					wypiszKolor(plik, naprezenie, naprezenieMax, glebia);
				}
				else
				{
					fprintf(plik, "%d %d %d", 0, 0, 0);
				}
				fprintf(plik, " ");
			}
			fprintf(plik, "\n");
		}
		fclose(plik);
		printf("\nZrobione!\n\n"); system("pause");
	}
	else
	{
		printf("\nB³¹d pliku!\n\n"); system("pause");
	}
}

void przekrojProstokatny()
{
	FILE* plik;
	double Ms, Mgx, Mgy, a, b, T, P, J, tau[2];
	double naprezenieMax = 0;
	double naprezenie;
	int x, y, wymiar, wymiarY, glebia;
	wczytywaniePodstawowychDanych(&Ms, &Mgx, &Mgy, &T, &P, &wymiar, &glebia);
	while (1)
	{
		printf("a[m]: ");	scanf("%lf", &a);
		if (a > 0) break;
		else printf("D³ugoœæ musi byæ dodatnia!\n");
	}
	while (1)
	{
		printf("b[m]: ");	scanf("%lf", &b);
		if (b > 0) break;
		else printf("D³ugoœæ musi byæ dodatnia!\n");
	}
	wymiarY = (int)(((double)wymiar)*b / a);
	tau[0] = 0; tau[1] = 0;
	if (a >= b)	J = K1(a / b)*a*b*b*b;
	else		J = K1(b / a)*b*a*a*a;
	for (y = 0; y < wymiarY; y++)
	{
		for (x = 0; x < wymiar; x++)
		{
			if(Ms!=0) tauObliczanie(&tau[0], a, b, Ms, J, (x - wymiar / 2)*a / wymiar, (wymiarY / 2 - y)*a / wymiar);
			naprezenie = 1/a/b*sqrt(pow(12 * (Mgx*(wymiarY / 2 - y)*a / wymiar /b/b - Mgy * (x - wymiar / 2)*a / wymiar/a/a) + P, 2) + 3 * (pow(a*b*tau[0],2) + pow(a*b*tau[1] + 6*T/b/b*(b*b/4 - (wymiarY / 2 - y)*(wymiarY / 2 - y)*a*a / wymiar / wymiar), 2)));
			if (naprezenie > naprezenieMax) naprezenieMax = naprezenie;
		}
	}
	plik = fopen("wizualizacja - przekrój prostok¹tny.ppm", "w");
	if (plik != NULL)
	{
		fprintf(plik, "P3\n%d %d\n%d\n", wymiar, wymiarY, glebia);
		for (y = 0; y < wymiarY; y++)
		{
			for (x = 0; x < wymiar; x++)
			{
			
				tauObliczanie(&tau[0], a, b, Ms, J, (x - wymiar / 2)*a / wymiar, (wymiarY / 2 - y)*a / wymiar);
				naprezenie = 1 / a / b * sqrt(pow(12 * (Mgx*(wymiarY / 2 - y)*a / wymiar / b / b - Mgy * (x - wymiar / 2)*a / wymiar / a / a) + P, 2) + 3 * (pow(a*b*tau[0], 2) + pow(a*b*tau[1] + 6 * T / b / b * (b*b / 4 - (wymiarY / 2 - y)*(wymiarY / 2 - y)*a*a / wymiar / wymiar), 2)));
				//naprezenie = 1 / a / b * sqrt(pow(12 * (Mgx*(wymiarY / 2 - y)*a / wymiar / b / b - Mgy * (x - wymiar / 2)*a / wymiar / a / a) + P, 2) + 3 * (/*ma byc a*b*tau poziome  pow(Ms*(x - wymiar / 2)*d / wymiar, 2)*/ +pow(/*tu ma byc a*b*tau pionowe -Ms * (wymiar / 2 - y)*d / wymiar*/ +6 * T / b / b * (b*b / 4 - (wymiarY / 2 - y)*(wymiarY / 2 - y)*a*a / wymiar / wymiar), 2)));
				wypiszKolor(plik, naprezenie, naprezenieMax, glebia);
				fprintf(plik, " ");
			}
			fprintf(plik, "\n");
		}
		fclose(plik);
		printf("\nZrobione!\n\n"); system("pause");
	}
	else
	{
		printf("\nB³¹d pliku!\n\n"); system("pause");
	}
}

int main()
{
	setlocale(LC_ALL, "polish_poland");
	przekrojProstokatny();
	
	return 0;
}
// zoptymalizowac przekrojKolowy()
// podzieliæ wczytywanie danych na sekcje tematyczne (wymiary geometryczne, obci¹¿enie, parametry grafiki)
// poprawic czarny pasek na dole
// dodac skrecanie dla wysokosci wiekszej od szerokosci