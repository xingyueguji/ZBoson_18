#include "./get_mc_eff.C"
#include "./get_all_bk_mc.C"
#include "./get_data.C"
#include "./get_ss_bk.C"
#include "./plot_normalize.C"

void run_all(bool doeff = 0, bool domc = 0, bool dodata = 0)
{

    cout << "Running all steps for preparation, this is eff " << doeff << " all mc background " << domc << " data " << dodata << endl;

    if (doeff)
        get_mc_eff();
    if (domc)
    {
        cout << "Doing mc 1" << endl;
        get_all_bk_mc(1);
        cout << "Doing mc 2" << endl;
        get_all_bk_mc(2);
        cout << "Doing mc 3" << endl;
        get_all_bk_mc(3);
    }
    if (dodata)
    {
        cout << "Doing data" << endl;
        get_data();
    }

    cout << "Doing same sign " << endl;

    get_ss_bk();

    plot_normalize(1);
    plot_normalize(2);
    plot_normalize(3);
    plot_normalize(4);
    // there's no opt == 5
    plot_normalize(6);
    plot_normalize(7);
    plot_normalize(8);

}