Shader "Unlit/ArTmo"
{
    Properties
    {
        _MainTex ("Texture", 2D) = "white" {}
        _ToneImg ("Tone Image", 2D) = "white" {}
        _ToneImgDs ("(unused)", 2D) = "white" {}
        _ToneImgLp ("LowPass Image", 2D) = "white" {}
        _TcsLut ("TCS LUT", 2D) = "white" {}
        _InDr ("In Dynamic Range", float) = 0.0
        _YDispPeak ("Y Peak Luminance", float) = 0.0
        _YComp ("Black Level", float) = 0.0
        _BeamsplitterComp ("Beamsplitter Comp", float) = 0.0

    }
    SubShader
    {
        Tags { "RenderType"="Opaque" }
        LOD 100

        Pass
        {
        GLSLPROGRAM // here begins the part in Unity's GLSL

        #version 130

        #extension GL_ARB_texture_rectangle : enable
        #extension GL_EXT_gpu_shader4 : enable


        uniform highp sampler2D _MainTex;
        uniform highp sampler2D _ToneImg;
        uniform highp sampler2D _TcsLut;  // tone, color-correction, sharpenning LUT2

        uniform highp float _InDr;        // Input dynamic range (of the tone curve) in the log2 units
        uniform highp float _YDispPeak;  // The peak luminance of the display
        uniform highp float _YComp;       // The black level we compensate for

        uniform highp float _BeamsplitterComp; // Compensation for the haploscope beamsplitter reflectivity

        uniform mat4x4 _M_native2display; // From native RGB to BT.2020 (required for the haploscope)

        uniform highp sampler2D _ToneImgLp; // low-pass image for unsharp masking


        // unusued
        uniform highp sampler2D _ToneImgDs; // downscaled (zones) tone image




         #ifdef VERTEX // here begins the vertex shader

         out vec4 texCoord;
         out vec4 texCoord1;
         out vec4 texCoord4;
         void main() // all vertex shaders define a main() function
         {
            gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
               // this line transforms the predefined attribute
               // gl_Vertex of type vec4 with the predefined
               // uniform gl_ModelViewProjectionMatrix of type mat4
               // and stores the result in the predefined output
               // variable gl_Position of type vec4.

            //texCoord.x = gl_MultiTexCoord0.x;
            texCoord = gl_MultiTexCoord0;
            texCoord1 = gl_MultiTexCoord0;
            texCoord4 = gl_MultiTexCoord0;
         }

         #endif // here ends the definition of the vertex shader


         #ifdef FRAGMENT // here begins the fragment shader


         in vec4 texCoord;
         in vec4 texCoord1;
         in vec4 texCoord4;        // sampler1D coordinates are aligned with the edges of the texels rather than their centers
        // This function corrects for thayt
        float correct_tex_coord( float x )
        {
            const float ts = 33.0; // How many texels does the texture have
            return x*(ts-1.)/ts + 0.5/ts;
        }

        // sampler3D coordinates are aligned with the edges of the texels rather than their centers
        // This function corrects for thayt
        vec3 correct_tex_coord_3d( vec3 x )
        {
            const float ts = 9.0; // How many texels does the texture have
            return x*(ts-1.)/ts + 0.5/ts;
        }


        void main()
        {

            float tone_in = texture(_ToneImg, texCoord1.st).x;
            float tone_norm = correct_tex_coord(tone_in/_InDr + 1.);
            vec3 tcs = textureLod(_TcsLut, vec2(0.0, 0.0), 0.).rgb; // (t)one, colour (c)orrection and (s)harpenning
            float tone_out = tcs.x;

            float de = tcs.z; // detail enhancement factor
            float base = texture(_ToneImgLp, texCoord4.st).x; // the low-pass of the tone mapped frame
            float detail = tone_out-base;
            tone_out = clamp( base + detail * de, -20., 0.);

            const float min_y = 1./8191.; // 1/(2^13-1) for 13-bit relative linear values
            vec3 rgb_in = (texture(_MainTex, texCoord.st).rgb-1.)*13.;  // log2(RGB)

            float s = tcs.y;  // Colour correction factor


            // Schlick with the exponent but in the log2 space
            // RGB_out is between 0 and 1, linear relative colour space
            vec3 rgb_out = (rgb_in-tone_in)*s+tone_out;

            vec3 RGB_out = min( (pow(vec3(2.), rgb_out) - min_y)/(1.-min_y), vec3(1.))*  mat3(_M_native2display);


            //mat3x3(vec3(0.753879233459602, 0.198671336179777, 0.0475975224268294),
            // vec3(0.0457648757737765, 0.941678347036237, 0.0125075575950154),
            // vec3(-0.00121109595995890, 0.0175939449477859, 0.983523842649768));

            gl_FragColor.rgb = max( RGB_out*(_YDispPeak+_YComp) - _YComp, vec3(0.005) ) * _BeamsplitterComp;
            // gl_FragColor.rgb = max( RGB_out, vec3(0.005) ) * _BeamsplitterComp;
            // gl_FragColor.rgb /= vec3(_YDispPeak);
            gl_FragColor.a = 1.;
        }


         #endif // here ends the definition of the fragment shader

         ENDGLSL // here ends the part in GLSL
        }
    }
}
