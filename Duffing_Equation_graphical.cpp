#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

// Duffing-Gleichung: y'' + δ*y' + α*y + β*y^3 = γ*cos(ω*t)

// Funktion, die das Differential der Duffing-Gleichung berechnet
void duffingEquationgraph(double t, double y[], double dydt[]) {
    double alpha = 1.0; // Konstante α
    double beta = -1.0; // Konstante β
    double delta = 0.2; // Dämpfungskonstante δ
    double gamma = 0.3; // Amplitudenkonstante γ
    double omega = 1.0; // Frequenzkonstante ω

    dydt[0] = y[1]; // Erste Ableitung von y ist die zweite Komponente des y-Vektors
    dydt[1] = gamma * cos(omega * t) - delta * y[1] - alpha * y[0] - beta * pow(y[0], 3); // Zweite Ableitung von y
}

// Runge-Kutta-Methode vierter Ordnung zur Integration der Differentialgleichungen
void rungeKuttagraph(double t0, double y0[], double t, double h) {
    double k1[2], k2[2], k3[2], k4[2];
    double yTemp[2];
    double dydt[2];

    ofstream output("duffing_plot.plt"); // Datei zur Ausgabe der Gnuplot-Skript öffnen

    double tNow = t0;
    double yNow[2];
    yNow[0] = y0[0];
    yNow[1] = y0[1];

    output << "set title 'Duffing Equation'\n";
    output << "set xlabel 'Time'\n";
    output << "set ylabel 'Position'\n";
    output << "plot '-' using 1:2 with lines title 'Position'\n";

    while (tNow < t) {
        // Berechne y für den nächsten Schritt
        duffingEquationgraph(tNow, yNow, dydt);
        k1[0] = h * dydt[0];
        k1[1] = h * dydt[1];

        yTemp[0] = yNow[0] + 0.5 * k1[0];
        yTemp[1] = yNow[1] + 0.5 * k1[1];
        duffingEquationgraph(tNow + 0.5 * h, yTemp, dydt);
        k2[0] = h * dydt[0];
        k2[1] = h * dydt[1];

        yTemp[0] = yNow[0] + 0.5 * k2[0];
        yTemp[1] = yNow[1] + 0.5 * k2[1];
        duffingEquationgraph(tNow + 0.5 * h, yTemp, dydt);
        k3[0] = h * dydt[0];
        k3[1] = h * dydt[1];

        yTemp[0] = yNow[0] + k3[0];
        yTemp[1] = yNow[1] + k3[1];
        duffingEquationgraph(tNow + h, yTemp, dydt);
        k4[0] = h * dydt[0];
        k4[1] = h * dydt[1];

        yNow[0] += (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6.0;
        yNow[1] += (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6.0;

        // Ausgabe in Datei
        output << tNow << "\t" << yNow[0] << endl;

        // Aktualisiere t
        tNow += h;
    }

    output << "e\n"; // Gnuplot "end" Befehl

    output.close(); // Datei schließen
}

int main() {
    double t0 = 0.0; // Anfangszeit
    double y0[2] = {0.1, 0.0}; // Anfangsbedingungen: y(0) = 0.1, y'(0) = 0
    double tFinal = 100.0; // Endzeit
    double h = 0.01; // Schrittweite

    rungeKuttagraph(t0, y0, tFinal, h); // Aufruf der Runge-Kutta-Methode

    cout << "Simulation abgeschlossen. Gnuplot-Skript wurde erstellt ('duffing_plot.plt')." << endl;

    return 0;
}
