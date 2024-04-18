using ScottPlot;
using System;
using CsvHelper;
using System.Globalization;
using CsvHelper.Configuration;
using AForge.Math;
public class point
{
    public double x { get; set;}
    public double y { get; set;}
}

class Program
{
    static void Main(string[] args)
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
            var complexes = records.Select(sample => new Complex((double)sample.x, 0.0)).ToArray();
            // 将complexes数组 length 补齐到2的幂次方
            int length = complexes.Length;
            int n = 1;
            double fs = 1e6;
            double eps = 1e-4;
            while (n < length)
            {
                n *= 2;
            }
            Complex[] newComplexes = new Complex[n];
            for (int i = 0; i < length; i++)
            {
                newComplexes[i] = complexes[i];
            }
            for (int i = length; i < n; i++)
            {
                newComplexes[i] = new Complex(0.0, 0.0);
            }
            FourierTransform.FFT(newComplexes, FourierTransform.Direction.Forward);
            var tmp = newComplexes.Select(sample => sample.Magnitude).ToArray();
            for(int i=0; i<tmp.Length; i++)
            {
                if (tmp[i] < eps) tmp[i] = 0;
            }
            
            var peaks = FindPeaks(tmp);
            
            int N = n;
            double res = 0;
            
            var idx = peaks.First();
            var basis = tmp[idx] * 2; 
            Console.WriteLine("idx:{0}",idx);
            Console.WriteLine("basis:{0}",basis);        

            for(int i=idx+1; i<=N/2; i++)
            {
                var magn = complexes[i].Magnitude*2;
                if (magn > eps && i % idx == 0)
                {
                    res += magn*magn;
                }
            }
            Console.WriteLine("res:{0}",res);
            var thd = Math.Sqrt(res)/basis;
            Console.WriteLine(thd);

            
            // plt();
            // useAforgeNet();

            
        }
    }

    static private void useFFTSharp()
    {
        
        var records = new double[16];
        for(int i=0; i<16; i++)
        {
            records[i] = 3*Math.Sin(2*Math.PI*i/16)+4*Math.Cos(2*Math.PI*i/32);
        }
        
        // Calculate the FFT as an array of complex numbers
        System.Numerics.Complex[] spectrum = FftSharp.FFT.Forward(records);

        // or get the magnitude (units²) or power (dB) as real numbers
        double[] magnitude = FftSharp.FFT.Magnitude(spectrum);

        Console.WriteLine("magnitude length:{0}",magnitude.Length);
        foreach (var item in magnitude)
        {
            Console.WriteLine(item);
        }
        
    }

    static private void useAforgeNet()
    {
        int N = 128;
        var records = new double[N];
        for(int i=0; i<N; i++)
        {
            records[i] = 3*Math.Sin(2*Math.PI*i/16)+1.5*Math.Sin(2*Math.PI*i/8)+0.9*Math.Sin(2*Math.PI*i/4);
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
            if (Math.Abs(complexes[i].Magnitude) > 1e-6) rec_x[i] = complexes[i].Magnitude;
            else rec_x[i] = 0.0;
        }
        ScottPlot.Plot plot = new ScottPlot.Plot();
        plot.Add.Scatter(x, rec_x);
        plot.SavePng("fft_sample.png", 800, 600);
        // use scottplot to plot the original signal
        ScottPlot.Plot plt = new ScottPlot.Plot();
        plt.Add.Scatter(x, records);
        plt.SavePng("origin_sample.png", 800, 600);

        
        var basis = complexes[1].Magnitude * 2;
        var idx = 1;
        for(int i=1; i<N; i++)
        {
            if (complexes[i].Magnitude > 1e-6)
            {
                basis = complexes[i].Magnitude * 2;
                idx = i;
                break;
            }
        }

        for(int i=idx+1; i<=N/2; i++)
        {
            var magn = complexes[i].Magnitude*2;
            if (magn > 1e-6)
            {
                res += magn*magn;
            }
        }
        var thd = Math.Sqrt(res)/basis;
        Console.WriteLine("THD:{0}",thd);

        for(int i=0; i<16; i++)
        {
            var real = complexes[i].Re;
            var image = complexes[i].Im;
            var magnitude = complexes[i].Magnitude;
            if (Math.Abs(complexes[i].Re) < 1e-6) real = 0;
            if (Math.Abs(complexes[i].Im) < 1e-6) image = 0;
            if (Math.Abs(complexes[i].Magnitude) < 1e-6) magnitude = 0;

            Console.WriteLine("实部:"+real + " " + "虚部:"+image + "幅值:"+magnitude);
        }

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
            var rec_x  = records.Select(sample => sample.x).ToArray();
            double[] x = new double[rec_x.Length];
            for(int i=0; i<rec_x.Length; i++)
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
