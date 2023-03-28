#include <zephyr/math/cfd/rotate.h>

namespace zephyr { namespace math {

using namespace geom;

void Rotate::to_local(Vector3d &velocity, const Vector3d &no) {
    Vector3d normal = no.normalized();

    // @formatter:off

    // матрица поворота
    Matrix3d R;
    if (no.z() == 0.0) {
        // Синус и косинус угла вращения
        double c = normal.x();
        double s = normal.y();

        R <<  c ,  s , 0.0,
             -s ,  c , 0.0,
             0.0, 0.0, 1.0;
    }
    else {
        // Ось вращения
        Vector3d n(0.0, normal.z(), -normal.y());
        n.normalize();

        // Синус и косинус угла вращения
        double c = normal.x();
        double s = sqrt(1.0 - c * c);

        // матрица поворота
        R <<  c          , -s * n.z()                    , s * n.y()                    ,
              s * n.z()  ,  c + (1 - c) * n.y() * n.y()  , (1 - c) * n.y() * n.z()      ,
             -s * n.y()  ,  (1 - c) * n.y() * n.z()      , c + (1 - c) * n.z() * n.z()  ;

    }
    // @formatter:on

    velocity = R * velocity;
}

void Rotate::to_global(Vector3d &velocity, const Vector3d &no) {
    Vector3d normal = no.normalized();


    // @formatter:off

    // матрица поворота
    Matrix3d R;
    if (no.z() == 0.0) {
        // Синус и косинус угла вращения
        double c = normal.x();
        double s = normal.y();

        R <<  c , -s , 0.0,
              s ,  c , 0.0,
             0.0, 0.0, 1.0;
    }
    else {
        // Ось вращения
        Vector3d n(0.0, -normal.z(), normal.y());
        n.normalize();

        // Синус и косинус угла вращения
        double c = normal.x();
        double s = sqrt(1.0 - c * c);

        // матрица поворота
        R <<  c          , -s * n.z()                    , s * n.y()                    ,
              s * n.z()  ,  c + (1 - c) * n.y() * n.y()  , (1 - c) * n.y() * n.z()      ,
             -s * n.y()  ,  (1 - c) * n.y() * n.z()      , c + (1 - c) * n.z() * n.z()  ;
    }
    // @formatter:on

    velocity = R * velocity;
}

static void to_local(geom::Matrix3d& stress, const geom::Vector3d& normal) {

}

static void to_global(geom::Matrix3d& stress, const geom::Vector3d& normal) {

}

}
}