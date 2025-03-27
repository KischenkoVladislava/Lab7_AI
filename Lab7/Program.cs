using System;
using System.Linq;
using System.Collections.Generic;

class Program
{
	static void Main(string[] args)
	{
		Console.WriteLine("Выберите задание:");
		Console.WriteLine("1 - Задание 1 (4x4 матрица, методы EM, RGMM, AN, \"линия\")");
		Console.WriteLine("2 - Задание 2 (6x6 матрица, поиск выбросов, корректировка, \"треугольник\")");
		Console.Write("Выбор: ");
		string choice = Console.ReadLine();

		if (choice == "1")
		{
			RunAssignment1();
		}
		else if (choice == "2")
		{
			RunAssignment2();
		}
		else
		{
			Console.WriteLine("Неверный выбор.");
		}

		Console.WriteLine("Нажмите любую клавишу для выхода...");
		Console.ReadKey();
	}

	//  Задание 1 (4x4)
	static void RunAssignment1()
	{
		// МПС 4x4 (вариант 4)
		double[,] D4 = new double[,]
		{
			{ 1.0,   1.0,   5.0,    0.2  },
			{ 1.0,   1.0,   3.0,    0.2  },
			{ 0.2,   0.333, 1.0,    0.1429 },
			{ 5.0,   5.0,   7.0,    1.0  }
		};

		Console.WriteLine("\n--- Задание 1 ---");
		Console.WriteLine("Матрица парных сравнений (4x4):");
		PrintMatrix(D4);

		// Метод EM (главное собственное число и вектор)
		(double lambdaMax, double[] emVector) = PowerMethodWithEigenvalue(D4, 1000, 1e-8);
		Normalize(emVector);

		Console.WriteLine("\nEM - Главное собственное число lambdaMax = {0:F4}", lambdaMax);
		Console.WriteLine("EM - Весовой вектор:");
		PrintVectorAsPercents(emVector);

		// Вычислим CI, CR
		int n = D4.GetLength(0);
		double CI = (lambdaMax - n) / (n - 1);
		// MRCI(n=4) 0.89
		double MRCI_for_n4 = 0.89;
		double CR = CI / MRCI_for_n4;
		Console.WriteLine("CI = {0:F4}, CR = {1:F4}", CI, CR);

		// RGMM
		double[] rgmm = GeometricMeanMethod(D4);
		Console.WriteLine("\nRGMM - Весовой вектор:");
		PrintVectorAsPercents(rgmm);

		// GCI
		double GCI = ComputeGCI(D4, rgmm);
		Console.WriteLine("GCI = {0:F4}, порог ~ 0.3526 (n=4)", GCI);

		// AN
		double[] an = AdditiveNormalizationMethod(D4);
		Console.WriteLine("\nAN - Весовой вектор:");
		PrintVectorAsPercents(an);

		// HCI/HCR
		double HCI = ComputeGCI(D4, an); 
		double randomHCI_4 = 1.9; // условно
		double HCR = HCI / randomHCI_4;
		Console.WriteLine("HCI = {0:F4}, HCR = {1:F4} (считаем порог 0.08..0.1)", HCI, HCR);

		// Метод "линия"
		Console.WriteLine("\nМетод \"линия\" - (эталон = 4)");
		double[] line4 = LineMethod(D4, 3);
		PrintVectorAsPercents(line4);

		double k_c_accurate = ComputeSpectralCoefficientFull(D4); // для матрицы 4×4
		Console.WriteLine("Спектральный коэффициент согласованности k_c: {0:F4}", k_c_accurate);
	}

	//  Задание 2 (6x6)
	static void RunAssignment2()
	{
		// МПС 6x6, варианта 4 
		double[,] D6 = new double[,]
		{
			{ 1.0,    2.0,    3.0,    6.0,    1.0,    2.0    }, 
            { 0.5,    1.0,    1.0,    2.0,    4.0,    5.0    },
            { 0.3333, 1.0,    1.0,    2.0,    0.3333, 1.0    }, 
            { 0.1667, 0.5,    0.5,    0.5,    1.0,    0.1667 }, 
            { 1.0,    0.25,   3.0,    1.0,    1.0,    1.0    }, 
            { 0.5,    0.2,    1.0,    6.0,    1.0,    1.0    }  
        };

		Console.WriteLine("\n--- Задание 2 ---");
		Console.WriteLine("Матрица парных сравнений (6x6):");
		PrintMatrix(D6);

		// Проверим исходный CR (методом EM)
		(double lambdaMax6, double[] em6) = PowerMethodWithEigenvalue(D6, 2000, 1e-10);
		Normalize(em6);
		int n = 6;
		double CI = (lambdaMax6 - n) / (n - 1);
		// Для n=6 средний индекс соглас. MRCI(n=6) 1.25 (по таблице Саати)
		double MRCI_for_n6 = 1.24;
		double CR6 = CI / MRCI_for_n6;
		Console.WriteLine("\nИсходная оценка CR (метод EM):");
		Console.WriteLine("lambdaMax = {0:F4}, CI = {1:F4}, CR = {2:F4}", lambdaMax6, CI, CR6);
		Console.WriteLine("(Порог CR 0.1 для n=6)");

		// Поиск выбросов
		Console.WriteLine("\nПоиск выбросов");
		// Метод укороченной МПС (CI)
		FindOutliersByShortMPS(D6);

		// Метод корреляций
		FindOutliersByCorrelation(D6);

		// Метод Xi-square (или лог.отклонения)
		FindOutliersByXiSquare(D6);

		// Автоматическая корректировка (WGMM, WAMM) с alpha = 0.5..0.9
		Console.WriteLine("\n[2] Автоматическая корректировка (WGMM, WAMM)");
		double[] alphaList = new double[] { 0.5, 0.6, 0.7, 0.8, 0.9 };
		foreach (double alpha in alphaList)
		{
			Console.WriteLine($"\n--- WGMM, alpha={alpha} ---");
			AutoCorrectMatrix(D6, alpha, method: "WGMM");

			Console.WriteLine($"\n--- WAMM, alpha={alpha} ---");
			AutoCorrectMatrix(D6, alpha, method: "WAMM");
		}

		Console.WriteLine("\nМетод \"треугольник\":");
		double[] triangleAccurate = TriangleMethodAccurate(D6);
		PrintVectorAsPercents(triangleAccurate);

		Console.WriteLine("\nЗадание 2 завершено.");
	}

	//  Вспомогательные методы

	static (double eigenvalue, double[] eigenvector) PowerMethodWithEigenvalue(double[,] matrix, int maxIterations, double tol)
	{
		int n = matrix.GetLength(0);
		double[] vec = new double[n];
		for (int i = 0; i < n; i++) vec[i] = 1.0;
		Normalize(vec);

		double[] newVec = new double[n];
		double eigenvalue = 0.0;
		for (int iter = 0; iter < maxIterations; iter++)
		{
			// newVec = D * vec
			for (int i = 0; i < n; i++)
			{
				double sum = 0.0;
				for (int j = 0; j < n; j++)
				{
					sum += matrix[i, j] * vec[j];
				}
				newVec[i] = sum;
			}
			double newEig = RayleighQuotient(vec, newVec);
			Normalize(newVec);
			if (Difference(vec, newVec) < tol)
			{
				eigenvalue = newEig;
				break;
			}
			Array.Copy(newVec, vec, n);
			eigenvalue = newEig;
		}
		return (eigenvalue, newVec);
	}
	static double RayleighQuotient(double[] v, double[] Dv)
	{
		double num = 0, den = 0;
		for (int i = 0; i < v.Length; i++)
		{
			num += v[i] * Dv[i];
			den += v[i] * v[i];
		}
		return (den == 0) ? 0 : (num / den);
	}

	// RGMM
	static double[] GeometricMeanMethod(double[,] D)
	{
		int n = D.GetLength(0);
		double[] gm = new double[n];
		for (int i = 0; i < n; i++)
		{
			double prod = 1.0;
			for (int j = 0; j < n; j++)
			{
				prod *= D[i, j];
			}
			gm[i] = Math.Pow(prod, 1.0 / n);
		}
		Normalize(gm);
		return gm;
	}

	// AN
	static double[] AdditiveNormalizationMethod(double[,] D)
	{
		int n = D.GetLength(0);
		double[] colSum = new double[n];
		for (int j = 0; j < n; j++)
		{
			double s = 0;
			for (int i = 0; i < n; i++)
				s += D[i, j];
			colSum[j] = s;
		}
		double[] v = new double[n];
		for (int j = 0; j < n; j++)
		{
			v[j] = 1.0 / colSum[j];
		}
		Normalize(v);
		return v;
	}

	// "Линия"
	static double[] LineMethod(double[,] D, int etalonIndex)
	{
		int n = D.GetLength(0);
		double[] w = new double[n];
		for (int i = 0; i < n; i++)
		{
			w[i] = (i == etalonIndex) ? 1.0 : D[i, etalonIndex];
		}
		Normalize(w);
		return w;
	}

	// GCI 
	static double ComputeGCI(double[,] D, double[] w)
	{
		int n = D.GetLength(0);
		double sum = 0;
		for (int i = 0; i < n; i++)
		{
			for (int j = i + 1; j < n; j++)
			{
				double lnDij = Math.Log(D[i, j]);
				double lnWij = Math.Log(w[i] / w[j]);
				double diff = lnDij - lnWij;
				sum += diff * diff;
			}
		}
		double factor = 2.0 / (n * (n - 1));
		return factor * sum;
	}

	// Спектральный коэффициент согласованности 
	static double ComputeSpectralCoefficientFull(double[,] D)
	{
		int n = D.GetLength(0);
		double lnN = Math.Log(n);
		// Получаем "линейные" векторы для каждого эталона
		List<double[]> lineVectors = new List<double[]>();
		for (int h = 0; h < n; h++)
		{
			double[] v = LineMethod(D, h); // LineMethod должна вернуть нормированный вектор весов
			lineVectors.Add(v);
		}

		List<double> kValues = new List<double>();
		foreach (var v in lineVectors)
		{
			double entropy = 0.0;
			for (int i = 0; i < v.Length; i++)
			{
				if (v[i] > 1e-12) // избегаем log(0)
					entropy -= v[i] * Math.Log(v[i]);
			}
			double normalizedEntropy = entropy / lnN; 
													  
			double kLocal = 1.0 - normalizedEntropy;
			kValues.Add(kLocal);
		}
		// Итоговый спектральный коэффициент - минимум среди локальных (наихудшая согласованность)
		double kSpectral = kValues.Min();
		return kSpectral;
	}

	// Поиск выбросов (Задание 2)
	// Метод укороченной МПС (CI)
	static void FindOutliersByShortMPS(double[,] D)
	{
		int n = D.GetLength(0);
		Console.WriteLine("\nМетод укороченной МПС (CI):");
		double minCI = double.MaxValue;
		int worstIndex = -1;

		for (int k = 0; k < n; k++)
		{
			double[,] Dshort = RemoveRowCol(D, k);
			(double lamShort, _) = PowerMethodWithEigenvalue(Dshort, 1000, 1e-8);
			int nShort = n - 1;
			double CIshort = (lamShort - nShort) / (nShort - 1);
			// Просто выводим
			Console.WriteLine($"  Исключая систему S{k + 1}: lambdaMax={lamShort:F4}, CI={CIshort:F4}");
			if (CIshort < minCI)
			{
				minCI = CIshort;
				worstIndex = k;
			}
		}
		Console.WriteLine($"=> Наименьший CI при исключении S{worstIndex + 1}, CI={minCI:F4}");
	}
	static double[,] RemoveRowCol(double[,] D, int index)
	{
		int n = D.GetLength(0);
		double[,] R = new double[n - 1, n - 1];
		int rr = 0;
		for (int i = 0; i < n; i++)
		{
			if (i == index) continue;
			int cc = 0;
			for (int j = 0; j < n; j++)
			{
				if (j == index) continue;
				R[rr, cc] = D[i, j];
				cc++;
			}
			rr++;
		}
		return R;
	}

	// Метод корреляций 
	static void FindOutliersByCorrelation(double[,] D)
	{
		int n = D.GetLength(0);
		Console.WriteLine("\nМетод корреляций:");

		// Средняя корреляция строк
		double[] rowCorr = new double[n];
		for (int i = 0; i < n; i++)
		{
			double sumCorr = 0; int count = 0;
			for (int k = 0; k < n; k++)
			{
				if (k == i) continue;
				double c = CorrelationRows(D, i, k);
				sumCorr += c; count++;
			}
			rowCorr[i] = (count > 0) ? sumCorr / count : 0.0;
		}
		// Средняя корреляция столбцов
		double[] colCorr = new double[n];
		for (int j = 0; j < n; j++)
		{
			double sumCorr = 0; int count = 0;
			for (int k = 0; k < n; k++)
			{
				if (k == j) continue;
				double c = CorrelationCols(D, j, k);
				sumCorr += c; count++;
			}
			colCorr[j] = (count > 0) ? sumCorr / count : 0.0;
		}

		// Ищем минимальные
		double minRow = double.MaxValue; int minRowIndex = -1;
		for (int i = 0; i < n; i++)
		{
			if (rowCorr[i] < minRow) { minRow = rowCorr[i]; minRowIndex = i; }
		}
		double minCol = double.MaxValue; int minColIndex = -1;
		for (int j = 0; j < n; j++)
		{
			if (colCorr[j] < minCol) { minCol = colCorr[j]; minColIndex = j; }
		}
		Console.WriteLine("  Средние корреляции строк:");
		for (int i = 0; i < n; i++)
		{
			Console.WriteLine($"    Строка {i + 1}: {rowCorr[i]:F3}");
		}
		Console.WriteLine("  Средние корреляции столбцов:");
		for (int j = 0; j < n; j++)
		{
			Console.WriteLine($"    Столбец {j + 1}: {colCorr[j]:F3}");
		}
		Console.WriteLine($"=> Минимальная корр. у строки {minRowIndex + 1}, столбца {minColIndex + 1}");
		Console.WriteLine("   => выброс: d[" + (minRowIndex + 1) + "," + (minColIndex + 1) + "]");
	}
	// Коэффициент корреляции Пирсона между строками i, k
	static double CorrelationRows(double[,] D, int i, int k)
	{
		int n = D.GetLength(0);
		// Соберём векторы row_i, row_k (длиной n), но исключим диагональ
		List<double> x = new List<double>();
		List<double> y = new List<double>();
		for (int j = 0; j < n; j++)
		{
			if (j == i || j == k) continue; // иногда исключают самих себя
			x.Add(D[i, j]);
			y.Add(D[k, j]);
		}
		return PearsonCorrelation(x, y);
	}
	// Коэффициента корреляции Пирсона между столбцами j, k
	static double CorrelationCols(double[,] D, int j, int k)
	{
		int n = D.GetLength(0);
		List<double> x = new List<double>();
		List<double> y = new List<double>();
		for (int i = 0; i < n; i++)
		{
			if (i == j || i == k) continue;
			x.Add(D[i, j]);
			y.Add(D[i, k]);
		}
		return PearsonCorrelation(x, y);
	}
	static double PearsonCorrelation(List<double> X, List<double> Y)
	{
		if (X.Count < 2) return 0.0;
		double mx = X.Average();
		double my = Y.Average();
		double num = 0, denx = 0, deny = 0;
		for (int i = 0; i < X.Count; i++)
		{
			double dx = X[i] - mx;
			double dy = Y[i] - my;
			num += dx * dy;
			denx += dx * dx;
			deny += dy * dy;
		}
		double den = Math.Sqrt(denx * deny);
		if (den < 1e-14) return 0.0;
		return num / den;
	}

	// Метод Xi-square 
	static void FindOutliersByXiSquare(double[,] D)
	{
		int n = D.GetLength(0);
		Console.WriteLine("\nМетод Xi-square (лог.отклонения):");
		// Сначала найдём "теоретические" веса (допустим, RGMM)
		double[] wRGMM = GeometricMeanMethod(D);
		List<double> tValues = new List<double>();
		Dictionary<(int, int), double> mapT = new Dictionary<(int, int), double>();
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				if (i == j) continue;
				double diff = Math.Abs(Math.Log(D[i, j]) - Math.Log(wRGMM[i] / wRGMM[j]));
				tValues.Add(diff);
				mapT[(i, j)] = diff;
			}
		}
		double meanT = tValues.Average();
		double stdT = Math.Sqrt(tValues.Average(x => (x - meanT) * (x - meanT)));
		double thresholdHigh = meanT + 2 * stdT; 

		Console.WriteLine($"  Среднее t = {meanT:F4}, std = {stdT:F4}, порог ~ {thresholdHigh:F4}");
		Console.WriteLine("  Выходящие за порог (потенциальные выбросы):");
		foreach (var kv in mapT)
		{
			if (kv.Value > thresholdHigh)
			{
				Console.WriteLine($"    d[{kv.Key.Item1 + 1},{kv.Key.Item2 + 1}] => t={kv.Value:F4}");
			}
		}
	}

	// Автоматическая корректировка (WGMM, WAMM)
	static void AutoCorrectMatrix(double[,] Dorig, double alpha, string method = "WGMM")
	{
		double[,] D = (double[,])Dorig.Clone();
		int n = D.GetLength(0);

		// Начальное вычисление CR
		(double lam, double[] w) = PowerMethodWithEigenvalue(D, 1000, 1e-8);
		Normalize(w);
		double CI = (lam - n) / (n - 1);
		double MRCI = 1.24; // для n=6
		double CR = CI / MRCI;
		Console.WriteLine($"k=0, CR={CR:F3}");

		int maxIter = 20; 
		for (int k = 1; k <= maxIter; k++)
		{
			if (CR <= 0.1)
			{
				Console.WriteLine($"=> Достигли CR={CR:F3} <= 0.1 за {k - 1} итераций");
				break;
			}
			// Обновим D(k+1)
			UpdateMatrix(D, w, alpha, method);
			// Пересчитаем CR
			(lam, w) = PowerMethodWithEigenvalue(D, 1000, 1e-8);
			Normalize(w);
			CI = (lam - n) / (n - 1);
			CR = CI / MRCI;
			Console.WriteLine($"k={k}, CR={CR:F3}");
			if (k == maxIter)
			{
				Console.WriteLine("Достигли maxIter, CR={0:F3}", CR);
			}
		}
	}
	static void UpdateMatrix(double[,] D, double[] w, double alpha, string method)
	{
		int n = D.GetLength(0);
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				if (i == j)
				{
					D[i, j] = 1.0;
				}
				else
				{
					double ratio = w[i] / w[j];
					if (method == "WGMM")
					{
						double dOld = D[i, j];
						double dNew = Math.Pow(dOld, alpha) * Math.Pow(ratio, 1 - alpha);
						D[i, j] = dNew;
					}
					else
					{
						double dOld = D[i, j];
						double dNew = alpha * dOld + (1 - alpha) * ratio;
						D[i, j] = dNew;
					}
				}
			}
		}
	}
	static double[,] Submatrix(double[,] D, List<int> indices)
	{
		int m = indices.Count;
		double[,] sub = new double[m, m];
		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < m; j++)
			{
				sub[i, j] = D[indices[i], indices[j]];
			}
		}
		return sub;
	}

	// Вычисление итоговых весов по методу "треугольник"
	static double[] TriangleMethodAccurate(double[,] D)
	{
		int n = D.GetLength(0);
		// Список всех "линейных" векторов (каждый альтернативный эталон)
		List<double[]> lineVectors = new List<double[]>();
		for (int i = 0; i < n; i++)
		{
			double[] wLine = LineMethod(D, i);
			lineVectors.Add(wLine);
		}

		// Инициализируем итоговый вектор - простое среднее всех "линейных" векторов
		double[] vFinal = new double[n];
		for (int i = 0; i < n; i++)
		{
			double sum = 0.0;
			foreach (var vec in lineVectors)
			{
				sum += vec[i];
			}
			vFinal[i] = sum / n;
		}
		Normalize(vFinal);

		// Итерационный процесс для уточнения итогового вектора
		double tol = 1e-6;
		int maxIter = 100;
		int iter = 0;
		double[] vOld = new double[n];
		while (iter < maxIter)
		{
			Array.Copy(vFinal, vOld, n);
			// Для каждого "линейного" вектора вычисляем расстояние до текущего итогового
			double[] weights = new double[lineVectors.Count];
			double weightSum = 0.0;
			for (int i = 0; i < lineVectors.Count; i++)
			{
				double d = EuclideanDistance(lineVectors[i], vFinal);
				// Если расстояние очень мало, задаём большое значение веса (чтобы этот вектор сильно влиял)
				double w = (d < 1e-8) ? 1e8 : 1.0 / (d);
				weights[i] = w;
				weightSum += w;
			}
			// Обновляем итоговый вектор как взвешенное среднее всех "линейных" векторов
			double[] vNew = new double[n];
			for (int j = 0; j < n; j++)
			{
				double sum = 0.0;
				for (int i = 0; i < lineVectors.Count; i++)
				{
					sum += weights[i] * lineVectors[i][j];
				}
				vNew[j] = sum / weightSum;
			}
			Normalize(vNew);
			// Проверяем сходимость 
			double diff = EuclideanDistance(vFinal, vNew);
			Array.Copy(vNew, vFinal, n);
			if (diff < tol)
				break;
			iter++;
		}
		Console.WriteLine($"\nМетод \"треугольник\" сошелся за {iter} итераций.");
		return vFinal;
	}

	static void Normalize(double[] v)
	{
		double s = v.Sum();
		if (Math.Abs(s) < 1e-15) return;
		for (int i = 0; i < v.Length; i++)
			v[i] /= s;
	}
	static double Difference(double[] a, double[] b)
	{
		double sum = 0;
		for (int i = 0; i < a.Length; i++)
		{
			sum += Math.Abs(a[i] - b[i]);
		}
		return sum;
	}
	static void PrintMatrix(double[,] M)
	{
		int rows = M.GetLength(0);
		int cols = M.GetLength(1);
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				Console.Write($"{M[i, j],8:F4} ");
			}
			Console.WriteLine();
		}
	}
	static void PrintVectorAsPercents(double[] w)
	{
		for (int i = 0; i < w.Length; i++)
		{
			Console.WriteLine($"   S{i + 1}: {w[i] * 100.0:F2}%");
		}
	}
	static double EuclideanDistance(double[] a, double[] b)
	{
		double sum = 0.0;
		for (int i = 0; i < a.Length; i++)
		{
			double d = a[i] - b[i];
			sum += d * d;
		}
		return Math.Sqrt(sum);
	}
}
