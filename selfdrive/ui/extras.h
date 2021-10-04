#include <time.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>

static void ui_draw_extras_limit_speed(UIState *s)
{
    const UIScene *scene = &s->scene;
    cereal::CarControl::SccSmoother::Reader scc_smoother = scene->car_control.getSccSmoother();
    int activeNDA = scc_smoother.getRoadLimitSpeedActive();
    int limit_speed = scc_smoother.getRoadLimitSpeed();
    int left_dist = scc_smoother.getRoadLimitSpeedLeftDist();

    if(activeNDA > 0)
    {
        int w = 185;
        int h = 60;
        int x = (s->fb_w + (bdr_s*2))/2 - w/2 - bdr_s - 844;
        int y = 238;

        const char* img = activeNDA == 1 ? "img_nda" : "img_hda";
        ui_draw_image(s, {x, y, w, h}, img, 1.f);
    }

    if(limit_speed > 10 && left_dist > 0)
    {
        int w = 160;
        int h = 160;
        int x = (bdr_s*2) + 1690;
        int y = 260;
        char str[32];

        nvgBeginPath(s->vg);
        nvgRoundedRect(s->vg, x, y, w, h, 210);
        nvgStrokeColor(s->vg, nvgRGBA(255, 0, 0, 200)); //red
        nvgStrokeWidth(s->vg, 33);
        nvgStroke(s->vg);

        NVGcolor fillColor = nvgRGBA(255, 255, 255, 200);//white
        nvgFillColor(s->vg, fillColor);
        nvgFill(s->vg);

        nvgFillColor(s->vg, nvgRGBA(0, 0, 0, 250));//black

        nvgFontSize(s->vg, 120);
        nvgFontFace(s->vg, "sans-bold");
        nvgTextAlign(s->vg, NVG_ALIGN_CENTER | NVG_ALIGN_MIDDLE);

        snprintf(str, sizeof(str), "%d", limit_speed);
        nvgText(s->vg, x+w/2, y+h/2, str, NULL);

        nvgFillColor(s->vg, nvgRGBA(255, 255, 255, 250));//white
        
        nvgFontSize(s->vg, 90);

        if(left_dist >= 1000)
            snprintf(str, sizeof(str), "%.1fkm", left_dist / 1000.f);
        else
            snprintf(str, sizeof(str), "%dm", left_dist);

        nvgText(s->vg, x+w/2, y+h + 70, str, NULL);
    }
    else
    {
        auto controls_state = (*s->sm)["controlsState"].getControlsState();
        int sccStockCamAct = (int)controls_state.getSccStockCamAct();
        int sccStockCamStatus = (int)controls_state.getSccStockCamStatus();

        if(sccStockCamAct == 2 && sccStockCamStatus == 2)
        {
            int w = 200;
            int h = 200;
            int x = (bdr_s*2) + 300;
            int y = 80;
            char str[32];

            nvgBeginPath(s->vg);
            nvgRoundedRect(s->vg, x, y, w, h, 210);
            nvgStrokeColor(s->vg, nvgRGBA(255, 0, 0, 200));
            nvgStrokeWidth(s->vg, 30);
            nvgStroke(s->vg);

            NVGcolor fillColor = nvgRGBA(0, 0, 0, 50);
            nvgFillColor(s->vg, fillColor);
            nvgFill(s->vg);

            nvgFillColor(s->vg, nvgRGBA(255, 255, 255, 250));

            nvgFontSize(s->vg, 140);
            nvgFontFace(s->vg, "sans-bold");
            nvgTextAlign(s->vg, NVG_ALIGN_CENTER | NVG_ALIGN_MIDDLE);

            nvgText(s->vg, x+w/2, y+h/2, "CAM", NULL);
        }
    }
}

static void ui_draw_extras(UIState *s)
{
    ui_draw_extras_limit_speed(s);
}
