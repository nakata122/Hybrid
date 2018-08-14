namespace smallpaint_fixed {
void render(int id, int size, int spp, double refr_index);
}
namespace smallpaint_painterly {
void render(int id, int size, int spp, double refr_index);
}
namespace smallpaint_bvh {
void render(int id, int size, int spp, double refr_index);
}
namespace smallpaint_pssmlt {
void render(int id, int size, int spp, double refr_index);
}
namespace smallpaint_smallmedia {
void render(int id, int size, int spp, double refr_index);
}
namespace smallpaint_ppm {
void render(int id, int size, int spp, double refr_index);
}
namespace smallpaint_vrl {
void render(int id, int size, int spp, int vrlps, int scene, int bounces, double sigma_a, double sigma_s, double g, int sampling, double refr_index, int intensity, bool mediumRadiance, bool surfaceRadiance);
}

