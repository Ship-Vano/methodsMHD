#pragma once

//#include "LinOp.h";
//#include "FileIO.h";
//
//template<typename DT, typename F>
//bool HLLDSolve(F func, DT t_0, DT T, DT tau, vector<DT> u_0, string filename)
//{
//	string path = "OutputData\\" + filename;
//	ofstream fpoints(path);
//	cout << "log[INFO]: Starting EulerSolve" << endl;
//	cout << "log[INFO]: Opening a file \"" << filename << "\" to write..." << endl;
//	if (fpoints.is_open())
//	{
//		DT t_i = t_0;
//		vector<DT> y_i = u_0;
//		/*fpoints << t_i << endl;
//		writeVectorToFile(fpoints, y_i);*/
//		int ind = 0;
//		vector<DT> y_ipp = u_0;
//		writeVectorToFile(fpoints, t_i, y_i);
//		while (abs(T - t_i) >= 1e-8)
//		{
//			y_ipp = y_i + tau * func(t_i, y_i);
//			//fpoints << t_i << endl;
//			y_i = y_ipp;
//			t_i += tau;
//			writeVectorToFile(fpoints, t_i, y_i);
//		}
//		fpoints.close();
//		return true;
//	}
//	else
//		cout << "log[ERROR]: Couldn't open or create a file" << endl;
//	return false;
//}