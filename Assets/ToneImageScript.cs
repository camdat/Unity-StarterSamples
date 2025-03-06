using System.Collections;
using System.Linq;

using System.Collections.Generic;
using static Unity.Mathematics.math;
using UnityEngine;

public class ToneImageScript : MonoBehaviour
{


    public int TextureLength =1024;

    public Texture2D texture;

    public float initialTone = 0;
    public bool inverse = false;


    public static double[] GetTcClahe(double[] l_in, double[] H, double dest_dr)
    {
        double delta = l_in[1] - l_in[0];
        H = H.Select(h => h * dest_dr).ToArray();
        double cont_enh = 1;
        double thr = cont_enh * delta;
        for (int kk = 0; kk < 10; kk++)
        {
            bool[] overMax = H.Select(h => h >= thr).ToArray();
            double e_b = H.Where((h, index) => overMax[index]).Sum(h => h - thr);
            if (e_b == 0)
                break;
            for (int i = 0; i < H.Length; i++)
            {
                if (overMax[i])
                    H[i] = thr;
            }
            double[] w = (double[])H.Clone();
            for (int i = 0; i < w.Length; i++)
            {
                if (overMax[i])
                    w[i] = 0;
            }
            double sum_w = w.Sum();
            double epsilon = 1e-6;
            if (sum_w < epsilon)
            {
                bool[] validBins = overMax.Select(b => !b).ToArray();
                if (validBins.All(b => !b))
                {
                    for (int i = 0; i < validBins.Length; i++)
                        validBins[i] = true;
                }
                double validBinsSum = validBins.Count(b => b);
                w = validBins.Select(b => b ? 1.0 / validBinsSum : 0).ToArray();
            }
            else
            {
                w = w.Select(weight => weight / sum_w).ToArray();
            }
            for (int i = 0; i < H.Length; i++)
            {
                H[i] += w[i] * e_b;
            }
        }
        double[] tc = new double[H.Length + 1];
        tc[0] = 0;
        for (int i = 1; i < tc.Length; i++)
        {
            tc[i] = tc[i - 1] + H[i - 1];
        }
        for (int i = 0; i < tc.Length; i++)
        {
            tc[i] -= dest_dr;
        }
        return tc;
    }
    public static double[] ComputeToneCurve(double Y_display_peak, double[] l_in, double[] H, double Y_amb, double Y_comp)
    {
        double dest_dr = log2((Y_display_peak + Y_amb) / Y_comp);
        return GetTcClahe(l_in, H, dest_dr);
    }

    public static double[] GetColorCorrectionLut(double[] l_in, double[] tc)
    {
        double delta = l_in[1] - l_in[0];
        // Calculate slope
        double[] slope = new double[tc.Length - 1];
        for (int i = 0; i < slope.Length; i++)
        {
            slope[i] = (tc[i + 1] - tc[i]) / delta;
        }
        // Calculate slope_l
        double[] slope_l = new double[l_in.Length - 1];
        for (int i = 0; i < slope_l.Length; i++)
        {
            slope_l[i] = l_in[i] + delta / 2;
        }
        // Convolution with a simple moving average
        int ksize = 5;
        double[] kernel = Enumerable.Repeat(1.0 / ksize, ksize).ToArray();
        slope = Convolve(slope, kernel);
        double k1 = 1.6774;
        double k2 = 0.9925;
        double[] g = l_in;
        double[] c = Interpolate(slope_l, slope, Clamp(g, slope_l.First(), slope_l.Last()));
        double[] cc_lut = c.Select(ci => max((1 + k1) * pow(ci, k2) / (1 + k1 * pow(ci, k2)), 0.4)).ToArray();
        return cc_lut;
    }
    private static double[] Convolve(double[] input, double[] kernel)
    {
        int ksize = kernel.Length;
        int halfKsize = ksize / 2;
        double[] result = new double[input.Length];
        for (int i = 0; i < input.Length; i++)
        {
            double sum = 0;
            for (int j = 0; j < ksize; j++)
            {
                int index = i + j - halfKsize;
                if (index >= 0 && index < input.Length)
                {
                    sum += input[index] * kernel[j];
                }
            }
            result[i] = sum;
        }
        return result;
    }


    private static double[] Interpolate(double[] x, double[] y, double[] xi)
    {
        double[] yi = new double[xi.Length];
        for (int i = 0; i < xi.Length; i++)
        {
            int j = System.Array.FindIndex(x, val => val > xi[i]) - 1;
            if (j < 0) j = 0;
            if (j >= x.Length - 1) j = x.Length - 2;



            double x0 = x[j], x1 = x[j + 1];
            double y0 = y[0], y1 = y[0]; // changed to fix crash
            yi[i] = y0 + (y1 - y0) * (xi[i] - x0) / (x1 - x0);
        }
        return yi;
    }
    private static double[] Clamp(double[] values, double min_v, double max_v)
    {
        return values.Select(v => max(min_v, min(max_v, v))).ToArray();
    }

    public static double[] GetUnsharpMaskLut(double[] l_in, double Y_display_peak, double[] tc, double Y_amb)
    {
        double delta = l_in[1] - l_in[0];
        // Calculate tc_slope
        double[] tc_slope = new double[tc.Length - 1];
        for (int i = 0; i < tc_slope.Length; i++)
        {
            tc_slope[i] = (tc[i + 1] - tc[i]) / delta;
        }
        // Append the last element to tc_slope
        tc_slope = tc_slope.Concat(new[] { tc_slope.Last() }).ToArray();
        // Calculate Y_out
        double[] Y_out = tc.Select(t => Y_display_peak * pow(2, t)).ToArray();
        // Calculate weber_boost
        double[] weber_boost = new double[Y_out.Length];
        for (int i = 0; i < Y_out.Length; i++)
        {
            weber_boost[i] = (Y_out[i] + Y_amb) / (Y_out[i] * tc_slope[i]);
        }
        // Calculate de_out_lut
        double[] de_out_lut = weber_boost.Select(wb => min(wb, 1.6)).ToArray();
        return de_out_lut;
    }




    // Start is called before the first frame update
    void Start()
    {

        // grab the current object and attach the tone curve
        /*
            pp = [-13 0]; % Minimum and maximum values that we handle (in the log2 domain)
            N; % The number of segments
            l_in; % position of tone-curve nodes as log2 lum
            frame = 1;  % Frame index
            Nf = 3; % filter size
            Y_display_peak;
            Y_display_black;
            Y_peak_src = 200;
            k_max;
            tf;
            P = [];
            a = 10^-3.7767;    % noise parameter a (dynamic component)
            b = 10.^-6.0776;   % noise parameter b (static component)

            ignore_zeros = true;


            Y_display_peak = 100 % Display max luminance
            Y_display_black = 1 % The black level of the display alone (without the ambient light)
            k_max = 1.6 % Maximum detail enhancement
            options.Y_white {isnumeric} = Y_display_peak % The white point of the content
            options.histogram {ischar} = "contrast-weighted" % "plain", "contrast-weighted", "zonal"


            vt.N = 32; %round((vt.pp(2)-vt.pp(1))/0.25);
            vt.l_in = linspace( vt.pp(1), vt.pp(2), vt.N+1 );
            vt.tf = TempFilterFIR(60); %TempFilterIIR( 60, 5 );

            vt.Y_display_peak = Y_display_peak;
            vt.Y_display_black = Y_display_black;
            vt.k_max = k_max;
            vt.options = options;

            H_prev = 0.0312
            Y_amb = 50
            Y_comp = 12.5
            Y_display_peak = 100
            dest_dr = 3.5850
            H = 0.111852


            delta = vt.l_in(2)-vt.l_in(1);
            slope = diff(tc)/delta;
            slope_l = vt.l_in(1:(end-1))+delta/2;

            ksize = 5;
            slope = conv( slope, ones(1,ksize)/ksize, 'same' );

            k1 = 1.6774;
            k2 = 0.9925;
            g = vt.l_in;
            c = interp1(slope_l, slope, clamp(g,slope_l(1),slope_l(end))); % contrast
            cc_lut = max((1+k1)*c.^k2 ./ (1+k1*c.^k2), 0.4);
        */


        double[] slope_l = new double[33] {
            -13.0000f, -12.5938f, -12.1875f, -11.7812f, -11.3750f, -10.9688f, -10.5625f, -10.1562f, -9.7500f, -9.3438f, -8.9375f,
            -8.5312f, -8.1250f, -7.7188f, -7.3125f, -6.9062f, -6.5000f, -6.0938f, -5.6875f, -5.2812f, -4.8750f, -4.4688f,
            -4.0625f, -3.6562f, -3.2500f, -2.8438f, -2.4375f, -2.0312f, -1.6250f, -1.2188f, -0.8125f, -0.4062f, 0.0000f
        };

        double[] tone_curve = ComputeToneCurve(100, slope_l, new double[] { 1/32, 1/32, 1/32, 1/32, 1/32, 1/32, 1/32, 1/32, 1/32, 1/32, 1/32, 1/32, 1/32, 1/32, 1/32, 1/32, 1/32, 1/32, 1/32, 1/32, 1/32, 1/32, 1/32, 1/32, 1/32, 1/32, 1/32, 1/32, 1/32, 1/32, 1/32, 1/32 }, 50, 12.5);

        // print tone_curve
        for (int i = 0; i < tone_curve.Length; i++) {
            Debug.Log("tone_curve[" + i + "] = " + tone_curve[i]);
        }

        double[] cc_lut = GetColorCorrectionLut(slope_l, tone_curve);

        // print cc_lut
        for (int i = 0; i < cc_lut.Length; i++) {
            Debug.Log("cc_lut[" + i + "] = " + cc_lut[i]);
        }

        double[] de_out_lut = GetUnsharpMaskLut(slope_l, 100, tone_curve, 50);

        // print de_out_lut
        for (int i = 0; i < de_out_lut.Length; i++) {
            Debug.Log("de_out_lut[" + i + "] = " + de_out_lut[i]);
        }

        double[] texData;

        // append tone_curve, cc_lut, and de_out_lut to texData
        texData = tone_curve.Concat(cc_lut).Concat(de_out_lut).ToArray();




        Renderer rend = GetComponent<Renderer>();

        // create a Texture2D that is the length of tonecurve, cc_lut, and de_out_lut, then set values in the R channel
        Texture2D texture1d = new Texture2D(tone_curve.Length + cc_lut.Length + de_out_lut.Length, 1, TextureFormat.RFloat, false);
        for (int i = 0; i < texData.Length; i++) {
            texture1d.SetPixel(i, 0, new Color((float)texData[i], 0,
            0, 1));
        }


        rend.material.SetTexture("_TcsLut", texture1d);


        // set the value of M_rgb2xyz
        /*   M = [0.636953507, 0.144619185, 0.168855854;
                        0.262698339, 0.678008766, 0.0592928953;
                        4.99407097e-17, 0.0280731358, 1.06082723]
        */
        Matrix4x4 M_rgb2xyz = new Matrix4x4( new Vector4(0.636953507f, 0.262698339f, 4.99407097e-17f, 0f),
                                                 new Vector4(0.144619185f, 0.678008766f, 0.0280731358f, 0f),
                                                 new Vector4(0.168855854f, 0.0592928953f, 1.06082723f, 0f),
                                                 new Vector4(0f, 0f, 0f, 0f));

        // set the value of M_p3toxyz
        /* M_p3_to_xyz = [0.4866,0.2657,0.1982;
                0.2290,0.6917,0.0793;
                0,0.0451,1.0437];
        */
        Matrix4x4 M_p3toxyz = new Matrix4x4( new Vector4(0.4866f, 0.2290f, 0f, 0f),
                                                 new Vector4(0.2657f, 0.6917f, 0f, 0f),
                                                 new Vector4(0.1982f, 0.0793f, 1.0437f, 0f),
                                                 new Vector4(0f, 0f, 0f, 0f));

        // inverse p3toxyz to get M_xyz2native
        Matrix4x4 M_xyz2native = M_p3toxyz.inverse;


        // set the value of M_xyz2display
        /*             M_xyz2display = [0.636953507, 0.144619185, 0.168855854;
                0.262698339, 0.678008766, 0.0592928953;
                4.99407097e-17, 0.0280731358, 1.06082723]
                */

        Matrix4x4 M_xyz2display = new Matrix4x4( new Vector4(0.636953507f, 0.262698339f, 4.99407097e-17f, 0f),
                                                 new Vector4(0.144619185f, 0.678008766f, 0.0280731358f, 0f),
                                                 new Vector4(0.168855854f, 0.0592928953f, 1.06082723f, 0f),
                                                 new Vector4(0f, 0f, 0f, 0f)).inverse;

        // set the value of M_native2xyz
        /*    M_native2xyz = [0.4866,0.2657,0.1982;
                0.2290,0.6917,0.0793;
                0,0.0451,1.0437];
                */

        Matrix4x4 M_native2xyz = new Matrix4x4( new Vector4(0.4866f, 0.2290f, 0f, 0f),
                                                 new Vector4(0.2657f, 0.6917f, 0.0451f, 0f),
                                                 new Vector4(0.1982f, 0.0793f, 1.0437f, 0f),
                                                 new Vector4(0f, 0f, 0f, 0f));

        Matrix4x4 M_native2display = M_xyz2display * M_native2xyz;


        rend.material.SetMatrix("_M_cont2native", M_xyz2native*M_rgb2xyz);


        RenderTexture buffer = new RenderTexture(
                               TextureLength,
                               TextureLength,
                               0,                            // No depth/stencil buffer
                               RenderTextureFormat.ARGB32   // Standard colour format
        );

        texture = new Texture2D(TextureLength,TextureLength,TextureFormat.ARGB32,true);

        MeshRenderer render = GetComponent<MeshRenderer>();
        Material material = render.sharedMaterial;

        Graphics.Blit(null, buffer, material);
        RenderTexture.active = buffer;           // If not using a scene camera

        texture.ReadPixels(
          new Rect(0, 0, TextureLength, TextureLength), // Capture the whole texture
          0, 0,                          // Write starting at the top-left texel
          false);                          // No mipmaps

        texture.Apply();


        // load the Tone Image into the ArTMO and Tone Extraction
        Renderer renderer = GameObject.Find("QuadTone").GetComponent<Renderer>();

        renderer.material.SetTexture("_ToneImg", texture);

        // now we do the same thing to QuadTone

        Texture2D texture_downsample = new Texture2D(TextureLength,TextureLength,TextureFormat.ARGB32,true);
        MeshRenderer render_downsample = GameObject.Find("QuadTone").GetComponent<MeshRenderer>();
        Material material_downsample = render_downsample.sharedMaterial;
        Graphics.Blit(null, buffer, material_downsample);
        RenderTexture.active = buffer;           // If not using a scene camera
        texture_downsample.ReadPixels(
          new Rect(0, 0, TextureLength, TextureLength), // Capture the whole texture
          0, 0,                          // Write starting at the top-left texel
          false);                          // No mipmaps
        texture_downsample.Apply();

        // now send them to upsampleedgestop
        Renderer renderer_upsample = GameObject.Find("QuadUpsample").GetComponent<Renderer>();
        renderer_upsample.material.SetTexture("_ToneImgDs", texture_downsample);
        renderer_upsample.material.SetTexture("_Image", texture);

        // lastly take edgestop and the others, and send them to ArTMO
        Texture2D texture_upsample = new Texture2D(TextureLength,TextureLength,TextureFormat.ARGB32,true);
        MeshRenderer render_upsample = GameObject.Find("QuadUpsample").GetComponent<MeshRenderer>();
        Material material_upsample = render_downsample.sharedMaterial;
        Graphics.Blit(null, buffer, material_upsample);
        RenderTexture.active = buffer;           // If not using a scene camera
        texture_upsample.ReadPixels(
          new Rect(0, 0, TextureLength, TextureLength), // Capture the whole texture
          0, 0,                          // Write starting at the top-left texel
          false);                          // No mipmaps
        texture_upsample.Apply();


        Renderer renderer_artmo = GameObject.Find("QuadArTMO").GetComponent<Renderer>();
        renderer_artmo.material.SetTexture("_ToneImg", texture);
        renderer_artmo.material.SetTexture("_ToneImgLp", texture_upsample);

        renderer_artmo.material.SetMatrix("_M_native2display", M_native2display);


        // cr rgb handled here
        // add M_rec709_rgbdisp
        /*    M_rec709_rgbdisp = [1, 0, 0;
        0, 1, 0;
        0, 0, 1];
        */

        Matrix4x4 M_rec709_rgbdisp = new Matrix4x4( new Vector4(1f, 0f, 0f, 0f),
                                                 new Vector4(0f, 1f, 0f, 0f),
                                                 new Vector4(0f, 0f, 1f, 0f),
                                                 new Vector4(0f, 0f, 0f, 0f));


        Renderer renderer_crrgb = GameObject.Find("QuadCrRGB").GetComponent<Renderer>();
        renderer_artmo.material.SetTexture("_ToneImg", texture); // ?, or texture from arTMO?
        renderer_artmo.material.SetMatrix("_M_rec709_rgbdisp", M_rec709_rgbdisp);



        // log success
        Debug.Log("Texture saved to ");


    }

    // Update is called once per frame
    void Update()
    {
      // add 1 to initialTone. If initialTone is greater than 1000, set inverse to true and subtract 1
      if (!inverse) {
        initialTone += 1;
      } else {
        initialTone -= 1;
      }
      if (initialTone >= 400 || initialTone <= 0) {
        inverse = !inverse;
      }

      Renderer renderer_artmo = GameObject.Find("QuadArTMO").GetComponent<Renderer>();
      renderer_artmo.material.SetFloat("_YDispPeak", initialTone);


    }
}
