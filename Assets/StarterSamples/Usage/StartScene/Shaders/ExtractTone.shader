Shader "Unlit/ExtractTone"
{
    Properties
    {
        _Image ("Texture", 2D) = "white" {}

        _TcsLut ("Tcs Lut", 2D) = "white" {}
        _InDr ("Input DR", Float) = 0.0

    }
    SubShader
    {
        Tags { "RenderType"="Opaque" }
        LOD 100

      Pass { // some shaders require multiple passes
         GLSLPROGRAM // here begins the part in Unity's GLSL
        #pragma multi_compile_local _ ADJUST_TO_LINEARSPACE



        #extension GL_ARB_texture_rectangle : enable
        #extension GL_EXT_gpu_shader4 : enable

        uniform sampler2D _Image;
        uniform sampler2D _TcsLut;  // tone, color-correction, sharpenning LUT

        uniform float _InDr;        // Input dynamic range (of the tone curve) in the log2 units

         #ifdef VERTEX // here begins the vertex shader

         out vec4 texCoord;
         void main() // all vertex shaders define a main() function
         {
            gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
               // this line transforms the predefined attribute
               // gl_Vertex of type vec4 with the predefined
               // uniform gl_ModelViewProjectionMatrix of type mat4
               // and stores the result in the predefined output
               // variable gl_Position of type vec4.

            texCoord = gl_MultiTexCoord0;

         }

         #endif // here ends the definition of the vertex shader

         #ifdef FRAGMENT // here begins the fragment shader


         in vec4 texCoord;
         float correct_tex_coord( float x )
        {

            const float ts = 33.0; // How many texels does the texture have

            return x*(ts-1.)/ts + 0.5/ts;

        }

        void main()
        {
            vec3 rgb_log = (texture(_Image, texCoord.st).rgb-1.)*13.;  // log2(RGB)


            float tone = max( max( rgb_log.r, rgb_log.g ), rgb_log.b );

            float tone_norm = correct_tex_coord(tone/_InDr + 1.); // To get values 0-1
            float tone_out = textureLod(_TcsLut, vec2(tone_norm, 0.), 0.).x;  // .x is the tone curve LUT

            gl_FragColor.rgb = vec3(tone, tone*tone, tone_out);
            gl_FragColor.a = 1.;
        }


         #endif // here ends the definition of the fragment shader

         ENDGLSL // here ends the part in GLSL
      }

    }
}
