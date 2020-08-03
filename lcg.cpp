#include <algorithm>
#include <iostream>
#include <map>
#include <math.h>
#include <vector>
#include <omp.h>
#include <ctime>

#include "lcg.h"

//Сложение по модулю
inline long long Generator::LCG::Add(const long long num1, const long long num2) {
  return (num1 + num2) % m;
}



//Вычитание по модулю
inline long long Generator::LCG::Subtr(const long long num1, const long long num2) {
  return (num1 - num2 + m) % m;
}



//Умножение по модулю
inline long long Generator::LCG::Mult(const long long num1, const long long num2) {
  return (num1 * num2) % m;
}



//Деление по модулю
inline long long Generator::LCG::Div(const long long num1, const long long num2) {
  return Mult(num1, Pow(num2, m - 2));
}



//Бинарное возведение в степень по модулю
inline long long Generator::LCG::Pow(const long long num, const long long pow) {
  if (pow == 0) {
    return 1;
  }
  if (pow & 1) {
    return Mult(Pow(num, pow - 1), num);
  }
  long long result = Pow(num, pow >> 1);
  return Mult(result, result);
}



//Преобразует значение X в eps, принадлежащее отрехку [0, 1]
inline double Generator::LCG::GetEps(const long long x) {
  return double(x) / m;
}



//Функция находит следующий после from элемент последовательности X
inline long long Generator::LCG::GetX(const long long from) {
  return Add(Mult(a, from), c);
}



//Функция находит следующий элемент последовательности eps
inline double Generator::LCG::Get() {
  x = GetX(x);
  return GetEps(x);
}



//Функция "перехода" X на k шагов от элемента from
long long Generator::LCG::GetXTransition(const long long from, const long long k) {
  long long a_to_k_pow = Pow(a, k);
  long long b = Subtr(a, 1);
  return Add(Mult(a_to_k_pow, from), Div(Mult(Subtr(a_to_k_pow, 1), c), b));
}



//Функция "перехода" eps на k шагов
double Generator::LCG::GetTransition(const long long k) {
  x = GetXTransition(x, k);
  return GetEps(x);
}



//Функция сравнения результатов процедуры "перехода" и последовательного
//вычисления элементов
bool Generator::LCG::Compare(const long long k) {
  double current_x = x;
  double sequential_result = 0;
  for (long long i = 0; i < k; ++i) {
    sequential_result = Get();
  }
  x = current_x;
  double result_by_transition = GetTransition(k);
  x = current_x;
  std::cout << "k = " << k << std::endl;
  std::cout << "Sequential result = " << sequential_result << std::endl;
  std::cout << "Result by transition = " << result_by_transition << std::endl;
  if (std::fabs(sequential_result - result_by_transition) < 1e-10) {
    std::cout << "The values are equal" << std::endl;
    std::cout << std::endl;
    return true;
  }
  std::cout << "The values are NOT equal!" << std::endl;
  std::cout << std::endl;
  return false;
}



double Generator::LCG::GoldsteinApproximation(const long long degrees_of_freedom) {
  std::vector<double> a{1.0000886, 0.4713941, 0.0001348028, -0.008553069,
    0.00312558, -0.0008426812, 0.00009780499};
  std::vector<double> b{-0.2237368, 0.02607083, 0.01128186, -0.01153761,
    0.005169654, 0.00253001, -0.001450117};
  std::vector<double> c{-0.01513904, -0.008986007, 0.02277679, -0.01323293,
    -0.006950356, 0.001060438, 0.001565326};
  //Уровень значимости
  double alpha = 0.90;
  double d = 2.0637 * pow(log(1 / (1 - alpha)) - 0.16, 0.4274) - 1.5774;
  double x2 = 0;
  for (int i = 0; i <= 6; ++i) {
    x2 += pow(degrees_of_freedom, -i / 2.0) * pow(d, i) * (a[i] + b[i]
      / degrees_of_freedom + c[i] / (degrees_of_freedom * degrees_of_freedom));
  }
  x2 = degrees_of_freedom * pow(x2, 3);
  return x2;
}



//Подсчет оценки критерия хи-квадрат
bool Generator::LCG::X2(const std::vector<long long>& culc_v, const std::vector<double>& theor_v,
  const long long degrees_of_freedom) {
  double x2 = 0;
  for (long long i = 0; i < culc_v.size(); ++i) {
    if (theor_v[i] == 0) {
      continue;
    }
    double num = culc_v[i] - theor_v[i];
    x2 += num * num / theor_v[i];
  }
  double theoretical_x2 = GoldsteinApproximation(degrees_of_freedom);
  std::cout << "Calculated value of X2 = " << x2 << std::endl;
  std::cout << "Theoretical value of X2 = " << theoretical_x2 << std::endl;
  if (x2 <= theoretical_x2 + 1e-10) {
    std::cout << "Test passed" << std::endl << std::endl;
    return true;
  }
  std::cout << "Test failed!" << std::endl << std::endl;
  x = x0;
  return false;
}



bool Generator::LCG::EquidistributionTest() {
  double start = clock();
  std::cout << "Equidistribution Test" << std::endl;
  x = x0;
  int d = 10;
  //v[i] - сколько раз встречается значение i в последовательности {y[i]},
  //где y[i] = eps[i] * d, eps - сгенерированные псевдослучайные числа
  //на отрезке [0, 1]
  std::vector<long long> v(d, 0);
  for (long long i = 0; i < n; ++i) {
    int y = Get() * d;
    ++v[y];
  }
  //p[i] - теоретическая вероятность появления значения i
  std::vector<double> p(d, 0);
  for (int i = 0; i < d; ++i) {
    p[i] = 1.0 / d;
    p[i] *= n;
  }
  x = x0;
  bool is_test_passed = X2(v, p, d - 1);
  double finish = clock();
//  std::cout << "Sequential runtime = " << (finish - start) / 1000.0
//    << std::endl;
//  std::cout << std::endl;
  return is_test_passed;
}



bool Generator::LCG::EquidistributionTestParallel() {
  double start = clock();
  std::cout << "Equidistribution Test Parallel" << std::endl;
  x = x0;
  int d = 10;
  //Число блоков последовательности, которые будут считаться параллельно
  //(= число нитей)
  const int num_threads = std::min(sqrt(n), 1e5);
  //v_threads[i][j] - сколько раз встречается значение j
  //в последовательности {y[j]} в i-ом блоке
  std::vector<std::vector<long long>> v_threads(num_threads,
    std::vector<long long>(d, 0));
  //Длина одного блока
  long long step = n / num_threads;
  if (n % num_threads) {
    ++step;
  }
  //from[i] - первый элемент i-го блока
  std::vector<long long> from(num_threads, 0);
  from[0] = GetX(x0);
  ++v_threads[0][GetEps(from[0]) * d];
  for (int i = 1; i < num_threads; ++i) {
    from[i] = GetXTransition(from[i - 1], step);
    ++v_threads[i][GetEps(from[i]) * d];
  }
#pragma omp parallel for num_threads(num_threads)
  for (int i = 0; i < num_threads; ++i) {
    for (int j = 1; j < step && i * step + j < n; ++j) {
      from[i] = GetX(from[i]);
      ++v_threads[i][GetEps(from[i]) * d];
    }
  }
  //Слияние результатов по блокам
  std::vector<long long> v(d, 0);
  for (int i = 0; i < d; ++i) {
    for (int j = 0; j < num_threads; ++j) {
      v[i] += v_threads[j][i];
    }
  }
  //p[i] - теоретическая вероятность появления значения i
  std::vector<double> p(d, 0);
  for (int i = 0; i < d; ++i) {
    p[i] = 1.0 / d;
    p[i] *= n;
  }
  x = x0;
  bool is_test_passed = X2(v, p, d - 1);
  double finish = clock();
//  std::cout << "Parallel runtime = " << (finish - start) / 1000.0 << std::endl;
//  std::cout << std::endl;
  return is_test_passed;
}



bool Generator::LCG::SerialTest() {
  double start = clock();
  std::cout << "Serial Test" << std::endl;
  x = x0;
  int d = 10;
  //Проверка критерия серий для пар элементов
  int k = 2;
  int size = d * d;
  long long num_groups = n / k;
  //v[i] - сколько раз встречается пара i
  std::vector<long long> v(size, 0);
  for (long long i = 0; i < num_groups; ++i) {
    int s = 0;
    for (int j = 0; j < k; ++j) {
      s = s * d + Get() * d;
    }
    ++v[s];
  }
  //p[i] - теоретическая вероятность появления пары i
  std::vector<double> p(size, 0);
  for (int i = 0; i < size; ++i) {
    p[i] = 1.0 / size;
    p[i] *= num_groups;
  }
  x = x0;
  bool is_test_passed = X2(v, p, size - 1);
  double finish = clock();
//  std::cout << "Sequential runtime = " << (finish - start) / 1000.0
//    << std::endl;
//  std::cout << std::endl;
  return is_test_passed;
}



bool Generator::LCG::SerialTestParallel() {
  double start = clock();
  std::cout << "Serial Test Parallel" << std::endl;
  x = x0;
  int d = 10;
  //Проверка критерия серий для пар элементов
  int k = 2;
  int size = d * d;
  long long num_groups = n / k;
  //Число блоков последовательности, которые будут считаться параллельно
  //(= число нитей)
  const int num_threads = std::min(sqrt(num_groups), 1e5);
  //v_threads[i][j] - сколько раз встречается пара j в i-ом блоке
  std::vector<std::vector<long long>> v_threads(num_threads,
    std::vector<long long>(size, 0));
  //Длина блока
  long long step = num_groups / num_threads;
  if (num_groups % num_threads) {
    ++step;
  }
  //from[i] - первый элемент i-го блока
  std::vector<long long> from(num_threads, 0);
  from[0] = GetX(x0);
  for (int i = 1; i < num_threads; ++i) {
    from[i] = GetXTransition(from[i - 1], step * k);
  }
#pragma omp parallel for num_threads(num_threads)
  for (int i = 0; i < num_threads; ++i) {
    for (int j = 0; j < step && i * step + j < num_groups; ++j) {
      int s = 0;
      for (int l = 0; l < k; ++l) {
        s = s * d + GetEps(from[i]) * d;
        from[i] = GetX(from[i]);
      }
      ++v_threads[i][s];
    }
  }
  //Слияниие результатов по блокам
  std::vector<long long> v(size, 0);
  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < num_threads; ++j) {
      v[i] += v_threads[j][i];
    }
  }
  //p[i] - теоретическая вероятность появления пары i
  std::vector<double> p(size, 0);
  for (int i = 0; i < size; ++i) {
    p[i] = 1.0 / size;
    p[i] *= num_groups;
  }
  x = x0;
  bool is_test_passed = X2(v, p, size - 1);
  double finish = clock();
//  std::cout << "Parallel runtime = " << (finish - start) / 1000.0 << std::endl;
//  std::cout << std::endl;
  return is_test_passed;
}



bool Generator::LCG::PokerTest() {
  double start = clock();
  std::cout << "Poker Test" << std::endl;
  x = x0;
  int d = 10;
  //Длина набора
  int length = 5;
  int categories = 7;
  long long num_groups = n / length;
  //v[i] - сколько наборов относится к i-ой категории, где
  //0-ая категория - все элементы различны
  //1-ая - 1 пара
  //2-ая - 2 пары
  //3-ья - 1 тройка
  //4-ая - 1 пара и 1 тройка
  //5-ая - 1 четверка
  //6-ая - все элементы одинаковые
  std::vector<long long> v(categories, 0);
  for (long long i = 0; i < num_groups; ++i) {
    std::map<int, int> y;
    for (int j = 0; j < length; ++j) {
      y[Get() * d]++;
    }
    if (y.size() == 5) {
      v[0]++;
      continue;
    }
    if (y.size() == 4) {
      v[1]++;
      continue;
    }
    if (y.size() == 1) {
      v[6]++;
      continue;
    }
    int max_count = 0;
    for (auto element : y) {
      max_count = std::max(max_count, element.second);
    }
    if (y.size() == 3) {
      if (max_count == 3) {
        v[3]++;
        continue;
      }
      v[2]++;
      continue;
    }
    if (max_count == 4) {
      v[5]++;
      continue;
    }
    v[4]++;
  }
  //p[i] - теоретическая вероятность появления i-го набора
  std::vector<double> p{0.30240, 0.50400, 0.10800, 0.07200, 0.00900, 0.00450,
    0.00010};
  for (int i = 0; i < p.size(); ++i) {
    p[i] *= num_groups;
  }
  x = x0;
  bool is_test_passed = X2(v, p, categories - 1);
  double finish = clock();
//  std::cout << "Sequential runtime = " << (finish - start) / 1000.0
//    << std::endl;
//  std::cout << std::endl;
  return is_test_passed;
}



bool Generator::LCG::PokerTestParallel() {
  double start = clock();
  std::cout << "Poker Test Parallel" << std::endl;
  x = x0;
  int d = 10;
  int length = 5;
  int categories = 7;
  int num_groups = n / length;
  //Число блоков последовательности, которые будут считаться параллельно
  //(= число нитей)
  const int num_threads = std::min(sqrt(num_groups), 1e5);
  //v_threads[i][j] - сколько наборов в i-ом блоке относится к категории j, где
  //0-ая категория - все элементы различны
  //1-ая - 1 пара
  //2-ая - 2 пары
  //3-ья - 1 тройка
  //4-ая - 1 пара и 1 тройка
  //5-ая - 1 четверка
  //6-ая - все элементы одинаковые
  std::vector<std::vector<long long>> v_threads(num_threads,
    std::vector<long long>(categories, 0));
  //Размер блока
  long long step = num_groups / num_threads;
  if (num_groups % num_threads) {
    ++step;
  }
  //from[i] - первый элемент i-го блока
  std::vector<long long> from(num_threads, 0);
  from[0] = GetX(x0);
  for (int i = 1; i < num_threads; ++i) {
    from[i] = GetXTransition(from[i - 1], step * length);
  }
#pragma omp parallel for num_threads(num_threads)
  for (int i = 0; i < num_threads; ++i) {
    for (int j = 0; j < step && i * step + j < num_groups; ++j) {
      std::map<int, int> y;
      for (int l = 0; l < length; ++l) {
        y[GetEps(from[i]) * d]++;
        from[i] = GetX(from[i]);
      }
      if (y.size() == 5) {
        v_threads[i][0]++;
        continue;
      }
      if (y.size() == 4) {
        v_threads[i][1]++;
        continue;
      }
      if (y.size() == 1) {
        v_threads[i][6]++;
        continue;
      }
      int max_count = 0;
      for (auto element : y) {
        max_count = std::max(max_count, element.second);
      }
      if (y.size() == 3) {
        if (max_count == 3) {
          v_threads[i][3]++;
          continue;
        }
        v_threads[i][2]++;
        continue;
      }
      if (max_count == 4) {
        v_threads[i][5]++;
        continue;
      }
      v_threads[i][4]++;
    }
  }
  //Слияние результатов по блокам
  std::vector<long long> v(categories, 0);
  for (int i = 0; i < categories; ++i) {
    for (int j = 0; j < num_threads; ++j) {
      v[i] += v_threads[j][i];
    }
  }
  //p[i] - теоретическая вероятность появления i-го набора
  std::vector<double> p{0.30240, 0.50400, 0.10800, 0.07200, 0.00900, 0.00450,
    0.00010};
  for (int i = 0; i < p.size(); ++i) {
    p[i] *= num_groups;
  }
  x = x0;
  bool is_test_passed = X2(v, p, categories - 1);
  double finish = clock();
//  std::cout << "Parallel runtime = " << (finish - start) / 1000.0 << std::endl;
//  std::cout << std::endl;
  return is_test_passed;
}



bool Generator::LCG::CouponCollectorsTest() {
  double start = clock();
  std::cout << "Coupon Collector's Test" << std::endl;
  x = x0;
  int d = 10;
  long long t = 10;
  //v[i], 0 <= i < t - 1 - число серий одинаковых элементов длины i + 1
  //v[t - 1] - число серий длины >= t
  std::vector<long long> v(t, 0);
  //Число рассмотренных интервалов (серий одинаковых элементов)
  long long num_intervals = 0;
  //previous_num - элементы из рассматриваемой серии
  //length - длина серии
  long long length = 1;
  for (long long i = 1, previous_num = Get() * d; i < n; ++i,
    ++length) {
    int y = Get() * d;
    if (y != previous_num) {
      ++v[std::min(length, t) - 1];
      ++num_intervals;
      length = 0;
      previous_num = y;
    }
  }
  ++v[std::min(length, t) - 1];
  ++num_intervals;
  //Вереятность того, что следующий элемент будет равен заданному 
  double p1 = 1.0 / d;
  //p[i] - теоретическая вероятность серии длины i
  std::vector<double> p(t, 0);
  //Степень числа p1
  double pow = 1;
  for (long long i = 1; i <= t; ++i, pow *= p1) {
    if (i < t) {
      p[i - 1] = pow * (1 - p1);
    }
    else {
      p[i - 1] = pow;
    }
    p[i - 1] *= num_intervals;
  }
  x = x0;
  bool is_test_passed = X2(v, p, t - 1);
  double finish = clock();
//  std::cout << "Sequential runtime = " << (finish - start) / 1000.0 << std::endl;
//  std::cout << std::endl;
  return is_test_passed;
}



bool Generator::LCG::CouponCollectorsTestParallel() {
  double start = clock();
  std::cout << "Coupon Collector's Test Parallel" << std::endl;
  x = x0;
  int d = 10;
  long long t = 10;
  //Число блоков последовательности, которые будут считаться параллельно
  //(= число нитей)
  const int num_threads = std::min(sqrt(n), 1e5);
  //v_threads[i][j], 0 <= j < t - 1 - число серий одинаковых элементов
  //длины j + 1 в i-ом блоке
  //v_threads[i][t - 1] - число серий длины >= t в i-ом блоке
  std::vector<std::vector<long long>> v_threads(num_threads,
    std::vector<long long>(t, 0));
  //Размер блока
  long long step = n / num_threads;
  if (n % num_threads) {
    ++step;
  }
  //from[i] - первый элемент блока
  std::vector<long long> from(num_threads, 0);
  from[0] = GetX(x0);
  for (int i = 1; i < num_threads; ++i) {
    from[i] = GetXTransition(from[i - 1], step);
  }
  //Крайние серии в каждом из блоков, где первый элемент пары - значение серии,
  //второй - длина серии
  std::vector<std::pair<int, long long>> first_intervals(num_threads, {-1, 0});
  std::vector<std::pair<int, long long>> last_intervals(num_threads, {-1, 0});
#pragma omp parallel for num_threads(num_threads)
  for (int i = 0; i < num_threads; ++i) {
    first_intervals[i].first = GetEps(from[i]) * d;
    int length = 0;
    while (int(GetEps(from[i]) * d) == first_intervals[i].first
      && length < step && i * step + length < n) {
      from[i] = GetX(from[i]);
      ++length;
    }
    first_intervals[i].second = length;
  }
  //Число рассмотренных интервалов (серий одинаковых элементов)
  long long num_intervals = 0;
#pragma omp parallel for reduction(+:num_intervals) num_threads(num_threads)
  for (int i = 0; i < num_threads; ++i) {
    //previous_num - элементы из рассматриваемой серии
    //length - длина серии
    long long previous_num = GetEps(from[i]) * d;
    from[i] = GetX(from[i]);
    long long length = 1;
    for (int j = first_intervals[i].second; j < step - 1 && i * step + j < n;
      ++j, ++length) {
      int y = GetEps(from[i]) * d;
      from[i] = GetX(from[i]);
      if (y != previous_num) {
        ++v_threads[i][std::min(length, t) - 1];
        ++num_intervals;
        length = 0;
        previous_num = y;
      }
    }
    if (i * step + first_intervals[i].second < n) {
      last_intervals[i] = std::make_pair(previous_num, length);
    }
  }
  //Слияние результатов по блокам за исключением граничных значений
  std::vector<long long> v(t, 0);
  for (int i = 0; i < t; ++i) {
    for (int j = 0; j < num_threads; ++j) {
      v[i] += v_threads[j][i];
    }
  }
  if (first_intervals[0].second == step) {
    //Весь первый блок состоит из однаковых элементов, необходимо сравнение
    //с началом второго
    last_intervals[0].first = first_intervals[0].first;
    last_intervals[0].second = step;
  }
  else {
    //Добавление серии из начала первого блока
    ++v[std::min(first_intervals[0].second, t) - 1];
    ++num_intervals;
  }
  for (int i = 0; i < num_threads - 1; ++i) {
    if (last_intervals[i].first == first_intervals[i + 1].first) {
      if (first_intervals[i + 1].second == step) {
        //Последняя серия i-го блока совпадает с первой серией i+1-го, причем
        //i+1 блок состоит из одинаковых элементов, необходимо сравнение с
        //началом i+2 блока
        last_intervals[i + 1].first = last_intervals[i].first;
        last_intervals[i + 1].second = last_intervals[i].second + step;
      }
      else {
        //Последняя серия i-го блока совпадает с первой серией i+1-го,
        //в вектор v добавляется результат слияние 2-х серий
        ++v[std::min(last_intervals[i].second + first_intervals[i + 1].second
          , t) - 1];
        ++num_intervals;
      }
    }
    else {
      //Последняя серия i-го блока не совпадает с первой серией i+1-го,
      //серия из i-го блока добавляется в вектор v
      ++v[std::min(last_intervals[i].second, t) - 1];
      ++num_intervals;
      if (first_intervals[i + 1].second == step) {
        //i+1 блок состоит из одинаковых элементов, необходимо сравнений
        //с началом i+2-го
        last_intervals[i + 1].first = first_intervals[i + 1].first;
        last_intervals[i + 1].second = step;
      }
      else {
        //Первыя серия i+1-го блока добавляется в вектор v
        ++v[std::min(first_intervals[i + 1].second, t) - 1];
        ++num_intervals;
      }
    }
  }
  //Добавление последней серии последнего блока, если она не пустая
  if (last_intervals[num_threads - 1].second != 0) {
    ++v[std::min(last_intervals[num_threads - 1].second, t) - 1];
    ++num_intervals;
  }
  //Вереятность того, что следующий элемент будет равен заданному 
  double p1 = 1.0 / d;
  //p[i] - теоретическая вероятность серии длины i
  std::vector<double> p(t, 0);
  //Степень числа p1
  double pow = 1;
  for (long long i = 1; i <= t; ++i, pow *= p1) {
    if (i < t) {
      p[i - 1] = pow * (1 - p1);
    }
    else {
      p[i - 1] = pow;
    }
    p[i - 1] *= num_intervals;
  }
  x = x0;
  bool is_test_passed = X2(v, p, t - 1);
  double finish = clock();
//  std::cout << "Sequential runtime = " << (finish - start) / 1000.0 << std::endl;
//  std::cout << std::endl;
  return is_test_passed;
}



bool Generator::LCG::Test() {
  if (EquidistributionTestParallel() && SerialTestParallel()
    && PokerTestParallel() && CouponCollectorsTestParallel()) {
    std::cout << "All tests passed!" << std::endl << std::endl;
    return true;
  }
  return false;
}



void Generator::LCG::ParallelTest() {
  EquidistributionTest();
  EquidistributionTestParallel();
  SerialTest();
  SerialTestParallel();
  PokerTest();
  PokerTestParallel();
  CouponCollectorsTest();
  CouponCollectorsTestParallel();
}
