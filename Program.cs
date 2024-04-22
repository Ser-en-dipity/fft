using ScottPlot;
using System;
using CsvHelper;
using System.Globalization;
using CsvHelper.Configuration;
using AForge.Math;
using static System.Math;

public class point
{
    public double x { get; set; }
    public double y { get; set; }
}

class Program
{
    static void Main(string[] args)
    {
        var config = new CsvConfiguration(CultureInfo.InvariantCulture)
        {
            HasHeaderRecord = false,
        };
        using (var reader = new StreamReader("C:/Users/小透明/Desktop/dq_err.csv"))
        using (var csv = new CsvReader(reader, config))
        {
            var records = csv.GetRecords<point>().ToList();
            var rec_x = records.Select(sample => sample.x).ToArray();
            int length = rec_x.Length;
            int n = 1;
            double eps = 1e-4;
            while (n < length)
            {
                n *= 2;
            }
            var new_rec_x = new double[n];
            for (int i = 0; i < length; i++)
            {
                new_rec_x[i] = rec_x[i];
            }
            for (int i = length; i < n; i++)
            {
                new_rec_x[i] = 0;
            }

            var window = new FftSharp.Windows.Hanning();
            double[] windowed = window.Apply(new_rec_x);

            System.Numerics.Complex[] spectrum = FftSharp.FFT.Forward(windowed);

            double[] magnitude = FftSharp.FFT.Magnitude(spectrum);

            // output the max magnitude
            var maxmag = 0.0;
            var idx = 0;
            for (int i = 0; i < magnitude.Length / 2; i++)
            {
                if (magnitude[i] > maxmag)
                {
                    idx = i;
                    maxmag = magnitude[i];
                }
            }

            Console.WriteLine("idx:{0}", idx);
            Console.WriteLine("basis:{0}", maxmag * 2);

            double[] x = new double[n];
            for (int i = 0; i < n; i++)
            {
                x[i] = i;
            }

            var magnitude_100 = magnitude.Take(500).ToArray();
            var x_100 = x.Take(500).ToArray();

            ScottPlot.Plot plot = new ScottPlot.Plot();
            plot.Add.Scatter(x_100, magnitude_100);
            plot.SavePng("fft.png", 800, 600);
            // use scottplot to plot the original signal
            ScottPlot.Plot plt = new ScottPlot.Plot();
            plt.Add.Scatter(x, rec_x);
            plt.SavePng("origin.png", 800, 600);

            var res = 0.0;
            for (int i = idx + 1; i <= n / 2; i++)
            {
                var magn = magnitude[i] * 2 * 2; // 汉宁窗幅值修正 and FFT 频谱对称
                if (magn > 1e-6 && i % idx == 0)
                {
                    res += magn * magn;
                }
            }
            var basis = maxmag * 2 * 2; // 汉宁窗幅值修正 and FFT 频谱对称
            var thd = Math.Sqrt(res) / basis;
            // console writeline thd
            Console.WriteLine("THD:{0}", thd);


            // plt();
            // useAforgeNet();
            // useFFTSharp();

        }
    }

    static private void useFFTSharp()
    {
        int N = 16384;
        var records = new double[N];
        for (int i = 0; i < N; i++)
        {
            records[i] = 8 * Sin(2 * PI * i / 1024) + 0.9 * Sin(2 * PI * i / 333);
        }
        var window = new FftSharp.Windows.Hanning();
        double[] windowed = window.Apply(records);

        // Calculate the FFT as an array of complex numbers
        System.Numerics.Complex[] spectrum = FftSharp.FFT.Forward(windowed);

        // or get the magnitude (units²) or power (dB) as real numbers
        double[] magnitude = FftSharp.FFT.Magnitude(spectrum);
        //  output the max magnitude
        var maxmag = 0.0;
        var idx = 0;
        for (int i = 0; i < magnitude.Length / 2; i++)
        {
            if (magnitude[i] > maxmag)
            {
                idx = i;
                maxmag = magnitude[i];
            }
        }
        Console.WriteLine("idx:{0}", idx);
        Console.WriteLine("basis:{0}", maxmag * 2);


        double[] x = new double[N];
        for (int i = 0; i < N; i++)
        {
            x[i] = i;
        }

        var magnitude_100 = magnitude.Take(100).ToArray();
        var x_100 = x.Take(100).ToArray();

        ScottPlot.Plot plot = new ScottPlot.Plot();
        plot.Add.Scatter(x_100, magnitude_100);
        plot.SavePng("fft_sample.png", 800, 600);
        // use scottplot to plot the original signal
        ScottPlot.Plot plt = new ScottPlot.Plot();
        plt.Add.Scatter(x, records);
        plt.SavePng("origin_sample.png", 800, 600);



    }

    static private void useAforgeNet()
    {
        int N = 16384;
        var records = new double[N];
        double xi = 0;
        while (xi < N / 2)
        {
            int index = (int)(2 * xi);
            records[index] = 9 * Sin(2 * PI * xi / 3000);
            xi += 0.5;
        }

        // create complex array
        var complexes = records.Select(sample => new Complex(sample, 0.0)).ToArray();
        FourierTransform.FFT(complexes, FourierTransform.Direction.Forward);
        double res = 0;
        // use scottplot to plot the fft result
        double[] rec_x = new double[N];
        double[] x = new double[N];
        for (int i = 0; i < N; i++)
        {
            x[i] = i;
            rec_x[i] = complexes[i].Magnitude;
        }

        var rec_100 = rec_x.Take(200).ToArray();
        var x_100 = x.Take(200).ToArray();
        ScottPlot.Plot plot = new ScottPlot.Plot();
        plot.Add.Scatter(x_100, rec_100);
        plot.SavePng("fft_sample.png", 800, 600);
        // use scottplot to plot the original signal
        ScottPlot.Plot plt = new ScottPlot.Plot();
        plt.Add.Scatter(x, records);
        plt.SavePng("origin_sample.png", 800, 600);


        var basis = 0.0;
        var idx = 0;
        for (int i = 1; i < N / 2; i++)
        {
            if (rec_x[i] > basis)
            {
                basis = rec_x[i];
                idx = i;
            }
        }
        basis = basis * 2;

        Console.WriteLine("idx:{0}", idx);
        Console.WriteLine("basis:{0}", basis);

        for (int i = idx + 1; i <= N / 2; i++)
        {
            var magn = complexes[i].Magnitude * 2;
            if (magn > 1e-6 && i % idx == 0)
            {
                res += magn * magn;
            }
        }
        var thd = Math.Sqrt(res) / basis;
        Console.WriteLine("THD:{0}", thd);

        // for (int i = 0; i < 16; i++)
        // {
        //     var real = complexes[i].Re;
        //     var image = complexes[i].Im;
        //     var magnitude = complexes[i].Magnitude;
        //     if (Math.Abs(complexes[i].Re) < 1e-6) real = 0;
        //     if (Math.Abs(complexes[i].Im) < 1e-6) image = 0;
        //     if (Math.Abs(complexes[i].Magnitude) < 1e-6) magnitude = 0;

        //     Console.WriteLine("实部:" + real + " " + "虚部:" + image + "幅值:" + magnitude);
        // }

    }

    static private void plt()
    {
        var config = new CsvConfiguration(CultureInfo.InvariantCulture)
        {
            HasHeaderRecord = false,
        };
        using (var reader = new StreamReader("C:/Users/i'CNC/Desktop/dq_err(1).csv"))
        using (var csv = new CsvReader(reader, config))
        {
            var records = csv.GetRecords<point>().ToList();
            // convert records to complex numbers
            var rec_x = records.Select(sample => sample.x).ToArray();
            double[] x = new double[rec_x.Length];
            for (int i = 0; i < rec_x.Length; i++)
            {
                x[i] = i;
            }

            ScottPlot.Plot plt = new ScottPlot.Plot();
            plt.Add.Scatter(x, rec_x);
            plt.SavePng("origin.png", 800, 600);
        }
    }
    static private void plot(double[] x, double[] y)
    {
        ScottPlot.Plot plt = new ScottPlot.Plot();
        plt.Add.Scatter(x, y);
        plt.SavePng("fft.png", 800, 600);
    }

    private static void FFT()
    {

    }

    public static int[] FindPeaks(double[] data)
    {
        double[] diff = new double[data.Length - 1];
        for (int i = 0; i < diff.Length; i++)
        {
            diff[i] = data[i + 1] - data[i];
        }
        int[] sign = new int[diff.Length];
        for (int i = 0; i < sign.Length; i++)
        {
            if (diff[i] > 0) sign[i] = 1;
            else if (diff[i] == 0) sign[i] = 0;
            else sign[i] = -1;
        }
        for (int i = sign.Length - 1; i >= 0; i--)
        {
            if (sign[i] == 0 && i == sign.Length - 1)
            {
                sign[i] = 1;
            }
            else if (sign[i] == 0)
            {
                if (sign[i + 1] >= 0)
                {
                    sign[i] = 1;
                }
                else
                {
                    sign[i] = -1;
                }
            }
        }
        List<int> result = new List<int>();
        for (int i = 0; i != sign.Length - 1; i++)
        {

            if (sign[i + 1] - sign[i] == -2)
            {
                result.Add(i + 1);
            }


        }
        return result.ToArray();//相当于原数组的下标
    }

}
