#include "DataHandler.h"
#include <cmath>

DataHandler::DataHandler(const char* filename) {
	fstream file(filename);
	if (!file) {
		cerr << "ERROR: could not open file '" << filename << "' for reading'" << endl;
		throw(-1);
	}
	// Geographic Information of Targets
	file >> num_targets;
	file >> depot1_loc.x >> depot1_loc.y;
	file >> depot2_loc.x >> depot2_loc.y;

	radii = new double[num_targets];
	for (int i = 0; i < num_targets; i++) {
		file >> radii[i];
	}
	target_locs = new Vertex[num_targets];
	for (int i = 0; i < num_targets; i++) {
		file >> target_locs[i].x >> target_locs[i].y;
	}
	// Task-related Information
	bdg_rewards = new double[num_targets];
	risk_thold = new double[num_targets];
	for (int i = 0; i < num_targets; i++) {
		file >> bdg_rewards[i];
	}
	for (int i = 0; i < num_targets; i++) {
		file >> risk_thold[i];
	}
}

DataHandler::~DataHandler() {
	delete radii;
	delete bdg_rewards;
	delete risk_thold;
	delete target_locs;
}

void DataHandler::print_data_info() {
	cout << "1. start depot:  " << depot1_loc.x << " " << depot1_loc.y << endl;
	cout << "2. end depot:  " << depot2_loc.x << " " << depot2_loc.y << endl;
	cout << "3. # of cirlces:  " << num_targets << endl;
	cout << "4. the lengths of radius are:  ";
	for (int i = 0; i < num_targets; i++)
		cout << radii[i] << "  ";

	cout << "\n";
	cout << "5. coordinates of targets' locations are:  \n";
	for (int i = 0; i < num_targets; i++)
		cout << target_locs[i].x << "\t" << target_locs[i].y << "\n";

	cout << "6. lower bounds of rewards:  ";
	for (int i = 0; i < num_targets; i++)
		cout << bdg_rewards[i] << "  ";

	cout << "\n";
	cout << "7. upper bounds of risks:  ";
	for (int i = 0; i < num_targets; i++)
		cout << risk_thold[i] << "  ";
}
