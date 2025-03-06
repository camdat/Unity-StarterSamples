Shader "Unlit/PrepareLogRGB"
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

        uniform mat4x4 _M_cont2native; // From content to native RGB




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

        vec3 srgb2lin( vec3 srgb )
        {
            const float t = 0.04045;
            const float a = 0.055;

            vec3 ltt = vec3(lessThanEqual(srgb, vec3(t)));
            vec3 rgb = srgb/12.92 * ltt + pow( (srgb+a)/(1.0+a), vec3(2.4) ) * (1.0-ltt);

            return rgb;
        }

        void main()
        {
            vec3 RGB_in_de = texture(_Image, texCoord.st).rgb;  // Display encoded
            vec3 RGB_in_cnt = srgb2lin( RGB_in_de ); // Linear (SDR content)
            vec3 RGB_in = clamp(mat3x3(_M_cont2native)*RGB_in_cnt, 0., 1.);  // From content to native RGB colour space

            const float min_y = 1./8191.; // 1/(2^13-1) for 13-bit relative linear values
            vec3 rgb_log = log2(RGB_in*(1.-min_y) + min_y);

            gl_FragColor.rgb = rgb_log / 13. + 1.;

            gl_FragColor.a = 1.;
        }


         #endif // here ends the definition of the fragment shader

         ENDGLSL // here ends the part in GLSL
      }

    Pass { // some shaders require multiple passes
         GLSLPROGRAM // here begins the part in Unity's GLSL


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
            //  vec3(pow(texture(_Image, texCoord.st).r, 2.0), pow(texture(_Image, texCoord.st).g, 2.0), pow(texture(_Image, texCoord.st).b, 2.0)); //
            vec3 rgb_log = (texture(_Image, texCoord.st).rgb-1.)*13.;  // log2(RGB)
            float tone =  max( max( rgb_log.r, rgb_log.g ), rgb_log.b );

            float tone_norm = correct_tex_coord(tone/_InDr + 1.); // To get values 0-1
            float tone_out = -5.1; //textureLod(_TcsLut, vec2(tone_norm, 0.), 0.).x;  // .x is the tone curve LUT

            gl_FragColor.rgb = vec3(tone, tone*tone, tone_out);
            gl_FragColor.a = 1.;
        }


         #endif // here ends the definition of the fragment shader

         ENDGLSL // here ends the part in GLSL
      }


    }
}
