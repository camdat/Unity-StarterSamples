Shader "Unlit/Pass"
{
    Properties
    {
        _Image ("Texture", 2D) = "black" {}

        _YDispPeak ("Lum Peak", Float) = 0.0
        _BeamsplitterComp ("Beamsplitter Comp", Float) = 0.0
        _FlipDirection ("Flip Direction", int) = 0

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
        uniform float _YDispPeak;  // The peak luminance of the display
        // uniform mat3 M_cont2display; // From content to native RGB
        uniform float _BeamsplitterComp; // Compensation for the haploscope beamsplitter reflectivity
        uniform highp int _FlipDirection; // Which way image should be flipped



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

            return vec3(rgb);
        }

        void main()
        {

            vec2 tex_cord = vec2( texCoord.s, texCoord.t ); // Non-modified texture coordinates
            if (_FlipDirection == 1) {
                float width = float(textureSize(_Image,0).x);
                tex_cord = vec2( width-texCoord.s, texCoord.t ); // Flip texture horizontally
            }
            if (_FlipDirection == 2) {
                float height = float(textureSize(_Image,0).y);
                tex_cord = vec2( texCoord.s, height-texCoord.t ); // Flip texture vertically
            }

            vec3 RGB_in_de = texture(_Image, tex_cord).rgb;  // Display encoded
            vec3 RGB_in_cnt = srgb2lin( RGB_in_de ); // Linear
            vec3 RGB_in = clamp( RGB_in_cnt, 0., 1.);  // From content to native RGB colour space

            gl_FragColor.rgb = RGB_in*_YDispPeak*_BeamsplitterComp;
            gl_FragColor.a = 1.;
        }


         #endif // here ends the definition of the fragment shader

         ENDGLSL // here ends the part in GLSL
      }

    }
}
