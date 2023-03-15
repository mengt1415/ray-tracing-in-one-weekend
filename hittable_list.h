#ifndef HITTABLE_LIST_H
#define HITTABLE_LIST_H
#include "rtweekend.h"
#include "hittable.h"
#include<vector>
#include<memory>

/*class hittable_list : public hittable {
public:
    hittable_list() {}
    hittable_list(shared_ptr<hittable> object) { add(object); }

    void clear() { objects.clear(); }
    void add(shared_ptr<hittable> object) { objects.push_back(object); }

    virtual bool hit(
        const ray& r, double t_min, double t_max, hit_record& rec) const override;

public:
    std::vector<shared_ptr<hittable>> objects;
};


bool hittable_list::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    hit_record temp_rec;
    auto hit_anything = false;
    auto closest_so_far = t_max;

    for (const auto& object : objects) {
        if (object->hit(r, t_min, closest_so_far, temp_rec)) {
            hit_anything = true;
            closest_so_far = temp_rec.t;
            rec = temp_rec;
        }
    }

    return hit_anything;
}*/
class hittable_list : public hittable {
public:
    hittable_list() {}
    hittable_list(hittable** l, int n) { list = l; list_size = n; }

    //void clear() { objects.clear(); }
    //void add(shared_ptr<hittable> object) { objects.push_back(object); }

    virtual bool hit(
        const ray& r, double t_min, double t_max, hit_record& rec) const;

    virtual bool bounding_box(double time0, double time1, aabb& output_box) const override;

public:
    hittable** list;
    //对于要生成的hittable__list，初始并不知道物体数组有多大
    //注意这里和my_sphere* slist的区别
    int list_size;
};


bool hittable_list::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    hit_record temp_rec;
    auto hit_anything = false;
    auto closest_so_far = t_max;

    for (int i = 0; i < list_size;i++) {
        if (list[i]->hit(r, t_min, closest_so_far, temp_rec)) {
            hit_anything = true;
            closest_so_far = temp_rec.t;
            rec = temp_rec;
        }
    }

    return hit_anything;
}
bool hittable_list::bounding_box(double time0, double time1, aabb& output_box) const 
{
    if (list_size < 1)
        return false;
    aabb temp_box;
    bool first_true = list[0]->bounding_box(time0, time1, temp_box);
    if (!first_true)
        return false;
    else
        output_box = temp_box;
    for (int i = 1; i < list_size; i++)
    {
        if (list[i]->bounding_box(time0, time1, temp_box))
        {
            output_box = surrounding_box(output_box, temp_box);
        }
        else
            return false;
    }
    return true;
}
#endif // !HITTABLE_LIST_H

