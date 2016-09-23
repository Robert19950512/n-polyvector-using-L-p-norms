#include <igl/avg_edge_length.h>
#include <igl/barycenter.h>
#include <igl/jet.h>
#include <igl/local_basis.h>
#include <igl/n_polyvector.h>
#include <igl/readDMAT.h>
#include <igl/readOBJ.h>
#include <igl/viewer/Viewer.h>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <igl/false_barycentric_subdivision.h>
#include <nanogui/formhelper.h>
#include <nanogui/screen.h>
#include <igl/streamlines.h>
#include "tutorial_shared_path.h"


// Input mesh
Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd VD;
Eigen::MatrixXi FD;
Eigen::MatrixXd Eout;
Eigen::MatrixXd K;
double p = 0.4;
double energy;
Eigen::VectorXi b(4);
double min = 0.0000001;
Eigen::MatrixXd bc(b.size(), 3);
Eigen::VectorXi samples;
std::stringstream te;
igl::StreamlineData sl_data;
igl::StreamlineState sl_state;




// Per face bases
Eigen::MatrixXd B1, B2, B3;
// Face barycenters
Eigen::MatrixXd B;

// Scale for visualizing the fields
double global_scale;

// Random length factor
double rand_factor = 5;

void representative_to_nrosy(
	const Eigen::MatrixXd &V,
	const Eigen::MatrixXi &F,
	const Eigen::MatrixXd &R,
	const int N,
	Eigen::MatrixXd &Y)
{
	using namespace Eigen;
	using namespace std;
	MatrixXd B1, B2, B3;

	igl::local_basis(V, F, B1, B2, B3);

	Y.resize(F.rows(), 3 * N);
	for (unsigned i = 0; i < F.rows(); ++i)
	{
		double x = R.row(i) * B1.row(i).transpose();
		double y = R.row(i) * B2.row(i).transpose();
		double angle = atan2(y, x);

		for (unsigned j = 0; j < N; ++j)
		{
			double anglej = angle + M_PI * double(j) / double(N);
			double xj = cos(anglej);
			double yj = sin(anglej);
			Y.block(i, j * 3, 1, 3) = xj * B1.row(i) + yj * B2.row(i);
		}
	}
}

Eigen::VectorXd random_constraints(const
	Eigen::VectorXd& b1, const
	Eigen::VectorXd& b2, int n)
{
	Eigen::VectorXd r(n * 3);
	for (unsigned i = 0; i<n; ++i)
	{
		double a = (double(rand()) / RAND_MAX) * 2 * M_PI;
		double s = 1 + ((double(rand()) / RAND_MAX)) * 5;
		Eigen::Vector3d t = s * (cos(a) * b1 + sin(a) * b2);
		r.block(i * 3, 0, 3, 1) = t;
	}
	return r;
}
bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifier)
{
	using namespace Eigen;
	using namespace std;
	if (key == 'A') {
		p=0.1;

	}
	else if (key == 'S') {
		p=0.5;
	}
	else if (key == 'D') {
		p = 0.9;
	}
	else if (key == 'F') {
		p = 2;
	}
	else {
		return false;
	}
	//std::cout << "p" << p << std::endl;
	viewer.data.clear();
	viewer.data.set_mesh(VD, FD);
	
	std::cout << "p" << p << std::endl;
	VectorXd k(FD.rows());
	k.setZero();
	K.resize(k.size(), 3);

	// for (unsigned i = 0; i<b.size(); ++i)
	//{



	// }
	
	Eigen::MatrixXd pvf;
	igl::n_polyvector(V, F, b, bc, pvf, Eout, p);
	Eigen::RowVector3d color = Eigen::RowVector3d::Ones();
	//std::cout << "pvf" << pvf << std::endl;
	igl::streamlines_init(V, F, pvf, true, sl_data, sl_state);
	for (int i = 0; i < 6; i++) {
		igl::streamlines_next(V, F, sl_data, sl_state);
		viewer.data.add_edges(sl_state.start_point, sl_state.end_point, color);
	}
	
	int Enum = Eout.rows();
	for (int eid = 0; eid < Enum; eid++) {
		int a = Eout(eid, 0);
		int b = Eout(eid, 1);

		int count = 0;
		for (int fid = 0; fid < FD.rows(); fid++) {
			for (int vid = 0; vid < FD.cols(); vid++) {
				if (FD(fid, vid) == a || FD(fid, vid) == b) {
					count++;

				}

			}


			if (count >= 2) {
				k[fid] = Eout(eid, 2);
			}
			count = 0;
		}
	}
	igl::jet(k, true, K);
	viewer.data.set_colors(K);
	//	cout << "k" << k << endl;
	for (int j = 0; j < k.size(); j++) {
		energy += k[j];
	}
	energy = energy / 2;
	cout << "total energy" << energy << endl;
	te << "p=" << p << "energy" << ": " << energy;
	Eigen::Vector3d m = V.colwise().maxCoeff();
	viewer.data.add_label(m, te.str());
	energy = 0;
	
//	Eigen::MatrixXd  temp_field2;
//	representative_to_nrosy(V, F, pvf, 3, temp_field2);
//	igl::streamlines_init(V, F, temp_field2, true, sl_data, sl_state);
	
	
	
	
	
	return false;
}






int main(int argc, char *argv[])
{
	using namespace Eigen;
	using namespace std;

	b << 4550, 2321, 5413, 5350;

	

	// Load a mesh in OBJ format
	igl::readOBJ(TUTORIAL_SHARED_PATH "/aircraft.obj", V, F);
	

	igl::false_barycentric_subdivision(V, F, VD, FD);
	igl::local_basis(V, F, B1, B2, B3);
	for (unsigned i = 0; i<b.size(); ++i)
	{
		VectorXd t = random_constraints(B1.row(b(i)), B2.row(b(i)), 1);
		bc.row(i) = t;
	}
	// Compute face barycenters
	igl::barycenter(V, F, B);

	// Compute scale for visualizing fields
	global_scale = .2*igl::avg_edge_length(V, F);
	

	igl::viewer::Viewer viewer;
	viewer.callback_key_down = &key_down;
	
	
	viewer.core.show_lines = false;
	viewer.launch();
	key_down(viewer, '2', 0);

}
